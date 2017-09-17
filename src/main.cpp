#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "uWS/uWS.h"
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "json.hpp"
#include "spline.h"
#include <map>

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
stringstream hasData(string s) {
	auto found_null = s.find("null");
	auto b1 = s.find_first_of("[");
	auto b2 = s.find_last_of("]");
	if (found_null != string::npos) {
		return stringstream();
	}
	else if (b1 != string::npos && b2 != string::npos) {
		stringstream tmp = stringstream();
		tmp.str(s.substr(b1, b2 - b1 + 1));
		return tmp;
	}
	return stringstream();
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for (int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x, y, map_x, map_y);
		if (dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y - y), (map_x - x));

	double angle = abs(theta - heading);

	if (angle > pi() / 4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

	int prev_wp;
	prev_wp = next_wp - 1;
	if (next_wp == 0)
	{
		prev_wp = maps_x.size() - 1;
	}

	double n_x = maps_x[next_wp] - maps_x[prev_wp];
	double n_y = maps_y[next_wp] - maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x + x_y*n_y) / (n_x*n_x + n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x, x_y, proj_x, proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000 - maps_x[prev_wp];
	double center_y = 2000 - maps_y[prev_wp];
	double centerToPos = distance(center_x, center_y, x_x, x_y);
	double centerToRef = distance(center_x, center_y, proj_x, proj_y);

	if (centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for (int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
	}

	frenet_s += distance(0, 0, proj_x, proj_y);

	return { frenet_s,frenet_d };

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp + 1) % maps_x.size();

	double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s - maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp] + seg_s*cos(heading);
	double seg_y = maps_y[prev_wp] + seg_s*sin(heading);

	double perp_heading = heading - pi() / 2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return { x,y };

}

vector<double> JMT(vector< double> start, vector <double> end, double T)
{
	/*
	Calculate the Jerk Minimizing Trajectory that connects the initial state
	to the final state in time T.

	INPUTS

	start - the vehicles start location given as a length three array
	corresponding to initial values of [s, s_dot, s_double_dot]

	end   - the desired end state for vehicle. Like "start" this is a
	length three array.

	T     - The duration, in seconds, over which this maneuver should occur.

	OUTPUT
	an array of length 6, each value corresponding to a coefficent in the polynomial
	s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

	EXAMPLE

	> JMT( [0, 10, 0], [10, 10, 0], 1)
	[0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
	*/

	Eigen::MatrixXd A = Eigen::MatrixXd(3, 3);
	A << T*T*T, T*T*T*T, T*T*T*T*T,
		3 * T*T, 4 * T*T*T, 5 * T*T*T*T,
		6 * T, 12 * T*T, 20 * T*T*T;

	Eigen::MatrixXd B = Eigen::MatrixXd(3, 1);
	B << end[0] - (start[0] + start[1] * T + .5*start[2] * T*T),
		end[1] - (start[1] + start[2] * T),
		end[2] - start[2];

	Eigen::MatrixXd Ai = A.inverse();

	Eigen::MatrixXd C = Ai*B;

	vector <double> result = { start[0], start[1], .5*start[2] };
	for (int i = 0; i < C.size(); i++)
	{
		result.push_back(C.data()[i]);
	}

	return result;

}

struct Vehicle {
	double x;
	double y;
	double d;
	double s;
	double v;
	int lane_number;
};

// Convert d-value to corresponding lane number
int convertDToLaneNumber(double d) {
	if (d < 4.0) { return 0; }
	if (d >= 4.0 and d < 8.0) { return 1; }
	if (d > 8.0) { return 2; }
	return -1;
}
// convert lane number to corresponding d-value
double convertLaneNumberToD(int lane_number) {
	if (lane_number == 0) { return 2.0; }
	if (lane_number == 1) { return 6.0; }
	if (lane_number == 2) { return 10.0; }
	return 14.0;
}

std::map<int, std::vector<Vehicle>> map_other_cars_to_their_lanes(const nlohmann::basic_json<> &sensor_fusion, int target_lane, double my_car_s)
{
	std::map<int, std::vector<Vehicle>> other_cars;
	for (auto sensor_values : sensor_fusion)
	{
		double vx = sensor_values[3];
		double vy = sensor_values[4];
		double v = sqrt(vx*vx + vy*vy);
		double s = sensor_values[5];
		float d = sensor_values[6];
		auto lane_number = convertDToLaneNumber(d);
		Vehicle other_car{
			sensor_values[1],
			sensor_values[2],
			d,
			s,
			v,
			lane_number
		};
		other_cars[other_car.lane_number].push_back(other_car);
	}
	return other_cars;
}

bool inline another_car_is_within_range(Vehicle my_car,
	int lane_number,
	double lower_bound_relative_to_my_car,
	double upper_bound_relative_to_my_car,
	const std::map<int,std::vector<Vehicle>> &other_cars)
{

	if (other_cars.find(lane_number) != other_cars.end())
	{
		auto absolute_upper = my_car.s + upper_bound_relative_to_my_car;
		auto absolute_lower = my_car.s + lower_bound_relative_to_my_car;
		for (const auto other_car : other_cars.at(lane_number))
		{
			if (other_car.s<absolute_upper&&other_car.s>absolute_lower)
			{
				return true;
			}
		}
	}
	return false;
}

bool inline can_change_to_left_lane(Vehicle my_car,
	const std::map<int, std::vector<Vehicle>> &other_cars,
	double relative_lower_bound_of_minimal_gap = 30,
	double relative_upper_bound_of_minimal_gap = -30)
{
	if (my_car.lane_number - 1 < 0)
	{
		return false;
	}
	else
	{
		return another_car_is_within_range(my_car,
			my_car.lane_number - 1,
			relative_lower_bound_of_minimal_gap,
			relative_upper_bound_of_minimal_gap,
			other_cars);
	};
}
bool inline can_change_to_right_lane(Vehicle my_car,
	const std::map<int, std::vector<Vehicle>> &other_cars,
	double relative_lower_bound_of_minimal_gap = 30,
	double relative_upper_bound_of_minimal_gap = -30)
{
	if (my_car.lane_number + 1 > 2)
	{
		return false;
	}
	else
	{
		return another_car_is_within_range(my_car,
			my_car.lane_number + 1,
			relative_lower_bound_of_minimal_gap,
			relative_upper_bound_of_minimal_gap,
			other_cars);
	};
}


int main() {
	uWS::Hub h;

	// Load up map values for waypoint's x,y,s and d normalized normal vectors
	vector<double> map_waypoints_x;
	vector<double> map_waypoints_y;
	vector<double> map_waypoints_s;
	vector<double> map_waypoints_dx;
	vector<double> map_waypoints_dy;

	// Waypoint map to read from
	string map_file_ = "../data/highway_map.csv";
	// The max s value before wrapping around the track back to 0
	double max_s = 6945.554;

	ifstream in_map_(map_file_.c_str(), ifstream::in);
	string line;


	while (getline(in_map_, line)) {
		istringstream iss(line);
		double x;
		double y;
		float s;
		float d_x;
		float d_y;
		iss >> x;
		iss >> y;
		iss >> s;
		iss >> d_x;
		iss >> d_y;
		map_waypoints_x.push_back(x);
		map_waypoints_y.push_back(y);
		map_waypoints_s.push_back(s);
		map_waypoints_dx.push_back(d_x);
		map_waypoints_dy.push_back(d_y);
	}

	//Start in lane 1
	int my_lane = 1;
	//The reference velocity
	double ref_vel = 0.0;

	h.onMessage([&ref_vel, &my_lane, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
		uWS::OpCode opCode) {
		// "42" at the start of the message means there's a websocket message event.
		// The 4 signifies a websocket message
		// The 2 signifies a websocket event
		//auto sdata = string(data).substr(0, length);
		//cout << sdata << endl;
		if (length && length > 2 && data[0] == '4' && data[1] == '2') {
			auto s = hasData(string(data));
			if (s.str() != "") {
				auto j = json::parse(s);
				string event = j[0].get<string>();
				if (event == "telemetry") {
//1: Get the current state of the main car 
					// j[1] is the data JSON object
					// Main car's localization Data
					double car_x = j[1]["x"];
					double car_y = j[1]["y"];
					double my_car_s = j[1]["s"];
					double car_d = j[1]["d"];
					double car_yaw = j[1]["yaw"];
					double car_speed = j[1]["speed"];
					Vehicle my_car{
						car_x,
						car_y,
						car_d,
						my_car_s,
						car_speed,
						my_lane
					};
//2: Get the previous path and the point up to which the car followed the previous path
					// Previous path data given to the Planner
					auto previous_path_x = j[1]["previous_path_x"];
					auto previous_path_y = j[1]["previous_path_y"];
					// Previous path's end s and d values 
					double end_path_s = j[1]["end_path_s"];
					double end_path_d = j[1]["end_path_d"];

//3: Get the sensor information representing the state of the other cars

					auto sensor_fusion = j[1]["sensor_fusion"];

//4: If there is a previous path, set the current s state of the main car to the last consumed point of the previous path, to guarantee a smooth transition

					int prev_size = previous_path_x.size();
					if (prev_size > 0)
					{
						my_car_s = end_path_s;
					}
					auto other_cars = map_other_cars_to_their_lanes(sensor_fusion, my_lane, my_car_s);
//5: Check if my lane is free and take actions if it is not
//5.1: If a car is in front, try to change lanes
					if (another_car_is_within_range(my_car,
						my_lane,
						0,
						30,
						other_cars))
					{
						auto new_target_lane = my_lane;
//5.1.1: Try to change to left lane first, if possible
						if (can_change_to_left_lane(my_car,
							other_cars))
						{
                            ///ToDo: Introduce a "prepare lane change" state
							new_target_lane = my_lane - 1;
						}
//5.1.2: If change to left lane is not possible, try to change to right lane
						else if (can_change_to_right_lane(my_car,
							other_cars))
						{
							///ToDo: Introduce a "prepare lane change" state
							new_target_lane = my_lane + 1;
						}
//5.2: Reduce speed to avoid crashing into car ahead
						ref_vel -= .224;
					}
//5.3: If lane is free, try to reach reference velocity
					else if (ref_vel < 49.5)
					{
						ref_vel += .224;
					}
//6: Create a spline as reference trajectory for the main car
//6.1: Create a list of widely spaced (x,y) waypoints, evenly spaced at a distance of 30m
					std::vector<double> ptsx;
					std::vector<double> ptsy;
//6.2: Create reference x, y, yaw states 
					double ref_x = car_x;
					double ref_y = car_y;
					double ref_yaw = deg2rad(car_yaw);

//6.3: If previous path is almost empty, use the car as starting reference
					if (prev_size < 2)
					{
						//Use two points that make the path tanget to the car
						double prev_car_x = car_x - cos(car_yaw);
						double prev_car_y = car_y - sin(car_yaw);
						ptsx.push_back(prev_car_x);
						ptsx.push_back(car_x);
						ptsy.push_back(prev_car_y);
						ptsy.push_back(car_y);
					}
//6.4: If the previous path is not almost empty, use the previous path's end point as starting reference
					else
					{
						//Redefine reference state as previous path and point
						ref_x = previous_path_x[prev_size - 1];
						ref_y = previous_path_y[prev_size - 1];
						double ref_x_prev = previous_path_x[prev_size - 2];
						double ref_y_prev = previous_path_y[prev_size - 2];
						ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);
						//Use two points t hat  make the path tangent to the previous path's and point
						ptsx.push_back(ref_x_prev);
						ptsx.push_back(ref_x);
						ptsy.push_back(ref_y_prev);
						ptsy.push_back(ref_y);
					}
//6.5: In Frenet add 3 evenly spaced waypoints 30 m apart from each others starting from the starting reference
					vector<double> next_wp0 = getXY(my_car_s + 30, (2 + 4 * my_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
					vector<double> next_wp1 = getXY(my_car_s + 60, (2 + 4 * my_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
					vector<double> next_wp2 = getXY(my_car_s + 90, (2 + 4 * my_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

					ptsx.push_back(next_wp0[0]);
					ptsx.push_back(next_wp1[0]);
					ptsx.push_back(next_wp2[0]);

					ptsy.push_back(next_wp0[1]);
					ptsy.push_back(next_wp1[1]);
					ptsy.push_back(next_wp2[1]);
//6.6: Convert points to car coordinates, to make the math easier
					for (int i = 0; i < ptsx.size(); i++)
					{
						double shift_x = ptsx[i] - ref_x;
						double shift_y = ptsy[i] - ref_y;
						ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y*sin(0 - ref_yaw));
						ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y*cos(0 - ref_yaw));
					}
//6.7: Create a spline
					tk::spline s;
					//set (x,y) points to the spline
					s.set_points(ptsx, ptsy);
//7: Calculate the actual (x,y) points along the spline, to pass to the planner
					vector<double> next_x_vals;
					vector<double> next_y_vals;
//7.1: Start with all the previous path points from last update
					for (int i = 0; i < previous_path_x.size(); i++)
					{
						next_x_vals.push_back(previous_path_x[i]);
						next_y_vals.push_back(previous_path_y[i]);
					}
//7.2: Calculate how to break up spline points so that we travel at our desired reference velocity
					double target_x = 30.0;
					double target_y = s(target_x);
					double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));
					double x_add_on = 0;
//7.3: Fill up the rest of our path planner after filling it with previous points, here we will always output 50 points
					for (int i = 1; i <= 50 - previous_path_x.size(); i++)
					{
						//Attention: Conversion to miles per hour
						double N = (target_dist / (.02*ref_vel / 2.24));
						double x_point = x_add_on + (target_x) / N;
						double y_point = s(x_point);
						x_add_on = x_point;
						double x_ref = x_point;
						double y_ref = y_point;
//7.4: Rotate back to global (map) coordinates
						x_point = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw));
						y_point = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw));
						x_point += ref_x;
						y_point += ref_y;
						next_x_vals.push_back(x_point);
						next_y_vals.push_back(y_point);
					}
//8: Pass the new path to the planner 
					json msgJson;
					msgJson["next_x"] = next_x_vals;
					msgJson["next_y"] = next_y_vals;
					auto msg = "42[\"control\"," + msgJson.dump() + "]";
					//this_thread::sleep_for(chrono::milliseconds(1000));
					ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
				}
			}
			else {
				// Manual driving
				std::string msg = "42[\"manual\",{}]";
				ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
			}
		}
	});

	// We don't need this since we're not using HTTP but if it's removed the
	// program
	// doesn't compile :-(
	h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
		size_t, size_t) {
		const std::string s = "<h1>Hello world!</h1>";
		if (req.getUrl().valueLength == 1) {
			res->end(s.data(), s.length());
		}
		else {
			// i guess this should be done more gracefully?
			res->end(nullptr, 0);
		}
	});

	h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
		std::cout << "Connected!!!" << std::endl;
	});

	h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
		char *message, size_t length) {
		ws.close();
		std::cout << "Disconnected" << std::endl;
	});

	int port = 4567;
	if (h.listen("0.0.0.0", port)) {
		std::cout << "Listening to port " << port << std::endl;
	}
	else {
		std::cerr << "Failed to listen to port" << std::endl;
		return -1;
	}
	h.run();
}
