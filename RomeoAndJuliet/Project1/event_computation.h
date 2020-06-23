#pragma once

#include<iostream>
#include<vector>
#include <algorithm>
#include "LOS.h"
#define INT_MAX 100000000
using namespace std;

LINE* minSumLine;
enum ROT {
	DEFAULT,
	CW,
	CCW
};

class EVENTS {
	int next_line_id;
	vector<vector<LINE*>> Queue;
	vector<int> shortest_path;
	vector<ROT> rotation; // true for CW, false for CCW
	SPT* spt[2]; //[0]: spt_s (s as root) , [1]: spt_t (t as root)
public:
	EVENTS() {}
	EVENTS(vector<int> _shortest_path, SPT* _spt_s, SPT* _spt_t)
	{
		next_line_id = 0;
		shortest_path = _shortest_path;
		for (int i = 0; i < shortest_path.size()-1; i++)
		{
			Queue.push_back(vector<LINE *> ());
		}
		spt[0] = _spt_s;
		spt[1] = _spt_t;
	}
	vector<int> get_shortest_path()
	{
		return shortest_path;
	}
	vector<vector<LINE*>> getQueue() {
		return Queue;
	}
	void sort_by_slope();
	void compute_path_events();
	void compute_boundary_events();
	void compute_bend_events();
	void compute_min_sum(void);
	double computeMinSum(void);
};

//sorts from small -> big
bool compare_slope(LINE* a, LINE* b)
{
	return (a->getSlope() < b->getSlope());
}



/* adds a LOS to the queue vector for every edge in the shortest path (s,t) */
void EVENTS::compute_path_events()
{
	for (int i = 0; i < shortest_path.size()-1; i++)
	{
		int prev = shortest_path[i];
		int cur = shortest_path[i + 1];
		
		PATH* line = new PATH(prev, cur,spt);
		Queue[i].push_back(line);
	}

	rotation.push_back(DEFAULT);
	for (int i = 1; i < shortest_path.size() - 1; i++)
	{
		double angle = calculate_angle_between(shortest_path[i], shortest_path[i + 1], shortest_path[i - 1], shortest_path[i]);
		//atan2 결과 양수이면 ccw 음수면 CW 0이면 전 결과 그대로
		if (angle > 0)
			rotation.push_back(CCW);
		else if (angle < 0)
			rotation.push_back(CW);
		else
			rotation.push_back(rotation.back());
	}
}

/* determines whether the line (*not vector) (CUR, P) is tangent to the path (PREV~CUR~NEXT) at vertex CUR */
bool is_tangent(int prev, int cur, int next, int p)
{
	float first = calculate_angle(prev, cur);
	float second = calculate_angle(cur, next);

	float p_cur = calculate_angle(p, cur);
	float cur_p = calculate_angle(cur, p);

	float angle1 = normalize_angle(first - p_cur);
	float angle2 = normalize_angle(p_cur - second);

	if (angle1*angle2 > 0)
		return true;
	
	angle1 = normalize_angle(first - cur_p);
	angle2 = normalize_angle(cur_p - second);

	if(angle1*angle2 > 0)
		return true;
	
	return false;
}

/* computes all boundary events
	sorts the boundary events by slope after inserting into the queue */
void EVENTS::compute_boundary_events()
{
	//search the tree for candidates 
	for (int i = 1; i < shortest_path.size() - 1; i++)
	{
		int prev = shortest_path[i - 1];
		int cur = shortest_path[i];
		int next = shortest_path[i + 1];

		//find the vertex in the spt and the direct children will be the candidate
		SPTnode* parent_s = spt[0]->get_node(cur);
		SPTnode* parent_t = spt[1]->get_node(cur);
		if (parent_s == NULL || parent_t ==NULL )
		{
			printf("couldn't find node in tree\n");
			exit(-1);
		}

		vector<SPTnode*> candidates(parent_s->get_children());
		int s_size = candidates.size();
		vector<SPTnode*> t_kids = parent_t->get_children();
		candidates.insert(candidates.end(), t_kids.begin(),t_kids.end());
		
		//then check the tangent thing...
		for (int j = 0; j < candidates.size(); j++)
		{
			int vertex_id = candidates[j]->get_id();
			if (vertex_id != next && vertex_id != shortest_path[i - 1])
			{
				if (is_tangent(prev,cur,next, vertex_id))
				{
					BOUNDARY* boundary = new BOUNDARY(cur, vertex_id,spt);
					Queue[i - 1].push_back(boundary);
				}
			}
		}
	}

	sort_by_slope();
}

void EVENTS::sort_by_slope()
{
	for (int i = 0; i < Queue.size()-1; i++)
	{
		double curS = Queue[i].front()->getSlope();
		double nextS = Queue[i + 1].front()->getSlope();

		std::sort(Queue[i].begin(), Queue[i].end(), compare_slope);
		vector<LINE*> newQueue;
		int curIdx = 0;
		for (; curIdx < Queue[i].size(); curIdx++)
		{
			if (curS == Queue[i][curIdx]->getSlope())
				break;
		}
		int nxtIdx = Queue[i].size();
		
		ROT rot = rotation[i+1];
		vector<LINE*>::iterator it = Queue[i].begin() + curIdx + 1;
		if (rot == CW) {
			newQueue.insert(newQueue.begin(), it, Queue[i].end());
			newQueue.insert(newQueue.end(), Queue[i].begin(), it);
			reverse(newQueue.begin(), newQueue.end());
			Queue[i] = newQueue;
		}
		else if (rot == CCW)
		{
			//newQueue.insert(newQueue.begin(), it-1, Queue[i].end());
			//newQueue.insert(newQueue.end(), Queue[i].begin(), it);
			//Queue[i] = newQueue;
			Queue[i].insert(Queue[i].end(), Queue[i].begin(), Queue[i].begin() + curIdx);
			Queue[i].erase(Queue[i].begin(), Queue[i].begin() + curIdx);
		}
	}
}

bool is_tangent_slope(double slope, double from, double to, ROT direction) {
	
	if (slope == from || slope == to)
		return false;
	bool inBetween = (direction == CW) == (from > to);
	if (inBetween)
		return (from - slope) * (to - slope) < 0;
	else
		return (from - slope) * (to - slope) > 0;
}
void EVENTS::compute_bend_events()
{
	int s_size = shortest_path.size();
	vector<int> prev = Queue[0][0]->getPath(0), cur, cur_t = Queue[0][0]->getPath(1), prev_t;
	double slope_prev = Queue[0][0]->getSlope(), slope_cur;
	for (int i = 0; i < Queue.size(); i++)
	{
		vector<LINE*> tempQueue;
		for (int j = 0; j < Queue[i].size(); j++)
		{
			cur = Queue[i][j]->getPath(0);
			prev_t = Queue[i][j]->getPath(1);
			slope_cur = Queue[i][j]->getSlope();

			//inspect paths with root _s
			int idx_prev = 0, idx_cur = 0;
			while (idx_cur != cur.size() || idx_prev != prev.size())
			{
				if (idx_cur != cur.size() && idx_prev != prev.size() && cur[idx_cur] == prev[idx_prev])
				{
					idx_cur++;
					idx_prev++;
				}
				else
				{
					BEND* bend;
					int rot = (j == 0) ? shortest_path[i] : shortest_path[i + 1];
					int idx = (j == 0) ? i : i + 1;
					if (idx_cur < cur.size()) {
						vector<int>::iterator it = find(prev.begin(), prev.end(), cur[idx_cur]);
						if (it == prev.end()) //type 1 (ii)
						{
							double slp = -1 / computeSlope(point_list[cur[idx_cur - 1]], point_list[cur[idx_cur]]);
							idx_cur++;
							if (is_tangent_slope(slp, slope_prev, slope_cur, rotation[idx]))
							{
								bend = new BEND(rot, cur[idx_cur - 2], cur[idx_cur - 1], spt, 0);
								if (bend->getType() != tERROR)
								{
									if (j == 0)
										Queue[i - 1].push_back(bend);
									else
										tempQueue.push_back(bend);
								}
							}
							/*
							bend = new BEND(rot, cur[idx_cur - 1], cur[idx_cur], spt,0);
							idx_cur++;
							if (bend->getType() != tERROR)
							{
								if (j == 0) {
									if (is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[i]))
										Queue[i - 1].push_back(bend);
								}
								else {
									if (is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[i + 1]))
										tempQueue.push_back(bend);
								}
							}*/
						}
					}
					if (idx_prev < prev.size()) {
						vector<int>::iterator it = find(cur.begin(), cur.end(), prev[idx_prev]);
						if (it == cur.end()) //type 2
						{
							double slp = -1 / computeSlope(point_list[prev[idx_prev - 1]], point_list[prev[idx_prev]]);
							idx_prev++;
							if (is_tangent_slope(slp, slope_prev, slope_cur, rotation[idx]))
							{
								bend = new BEND(rot, prev[idx_prev - 2], prev[idx_prev - 1], spt, 0);
								if (bend->getType() != tERROR)
								{
									if (j == 0)
										Queue[i - 1].push_back(bend);
									else
										tempQueue.push_back(bend);
								}
							}

							/*
							bend = new BEND(rot, prev[idx_prev - 1], prev[idx_prev], spt,0);
							idx_prev++;
							if (bend->getType() != tERROR)
							{
								if (j == 0) {
									if (is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[i]))
										Queue[i - 1].push_back(bend);
								}
								else {
									if (is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[i + 1]))
										tempQueue.push_back(bend);
								}
							}*/
						}
					}

				}
			}
			//bend events in terms of paths with root _t
			idx_prev = 0, idx_cur = 0;
			while (idx_cur != cur_t.size() || idx_prev != prev_t.size())
			{
				if (idx_cur != cur_t.size() && idx_prev != prev_t.size() && cur_t[idx_cur] == prev_t[idx_prev])
				{
					idx_cur++;
					idx_prev++;
				}
				else
				{
					BEND* bend;
					int rot = (j == 0) ? shortest_path[i] : shortest_path[i + 1];
					int idx = (j == 0) ? i : i + 1;

					if (idx_cur < cur_t.size()) {
						vector<int>::iterator it = find(prev_t.begin(), prev_t.end(), cur_t[idx_cur]);
						if (it == prev_t.end()) //type 1 (ii)
						{
							double slp = -1 / computeSlope(point_list[cur_t[idx_cur - 1]], point_list[cur_t[idx_cur]]);
							idx_cur++;
							if (is_tangent_slope(slp, slope_prev, slope_cur, rotation[idx])) {
								bend = new BEND(rot, cur_t[idx_cur-2], cur_t[idx_cur-1], spt,1);
								
								if (bend->getType() != tERROR) {
									if (j == 0)
										Queue[i - 1].push_back(bend);
									else
										tempQueue.push_back(bend);
								}
							}
							/*
							bend = new BEND(rot, cur_t[idx_cur - 1], cur_t[idx_cur], spt);
							idx_cur++;
							if (bend->getType() != tERROR)
							{
								if (j == 0) {
									if (is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[i]))
										Queue[i - 1].push_back(bend);
								}
								else {
									if (is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[i + 1]))
										tempQueue.push_back(bend);
								}
							}*/
						}
					}
					if (idx_prev < prev_t.size()) {
						vector<int>::iterator it = find(cur_t.begin(), cur_t.end(), prev_t[idx_prev]);
						if (it == cur_t.end()) //type 2
						{
							double slp = -1 / computeSlope(point_list[prev_t[idx_prev - 1]], point_list[prev_t[idx_prev]]);
							idx_prev++; 
							if (is_tangent_slope(slp, slope_prev, slope_cur, rotation[idx]))
							{
								bend = new BEND(rot, prev_t[idx_prev - 2], prev_t[idx_prev-1], spt,1);
								
								if (bend->getType() != tERROR)
								{
									if (j == 0)
										Queue[i - 1].push_back(bend);
									else
										tempQueue.push_back(bend);
								}
							}

							/*
							bend = new BEND(rot, prev_t[idx_prev - 1], prev_t[idx_prev], spt);
							idx_prev++;
							if (bend->getType() != tERROR)
							{
								if (j == 0) {
									if (is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[i]))
										Queue[i - 1].push_back(bend);
								}
								else {
									if (is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[i + 1]))
										tempQueue.push_back(bend);
								}
							}*/
						}
					}
				}
			}

			prev = cur;
			cur_t = prev_t;
			slope_prev = slope_cur;
		}
		Queue[i].insert(Queue[i].end(), tempQueue.begin(), tempQueue.end());
	}

	sort_by_slope();
}

double dist(double slope, int u, int v)
{
	Point U = point_list[u];
	Point V = point_list[v];

	double distance = abs(slope * U.get_x() - U.get_y() - slope * V.get_x() + V.get_y());
	distance /= sqrt(slope * slope + 1);
	return distance;
}

double getSlopeMinSum(double bound1,double bound2, ROT dir, int v, int u, int u_, double* sum)
{
	double minSlope = bound1;
	/* need to check three possible cases
	(i) boundary slopes
	(ii) absolute value turn slopes
	(iii) local optimum/maximum slopes */
	
	if (u == -1 || u_ == -1 || v == -1)
		printf("whatthe");
	double minSumDist = dist(bound1, u, v) + dist(bound1, u_, v);
	
	double distance = dist(bound2, u, v) + dist(bound2, u_, v);
	if (minSumDist > distance)
	{
		minSlope = bound2;
		minSumDist = distance;
	}

	//(ii) case
	Point U = point_list[u];
	Point U_ = point_list[u_];
	Point V = point_list[v];
	double c1 = U.get_x() - V.get_x();
	double c2 = V.get_y() - U.get_y();
	double c_1 = U_.get_x() - V.get_x();
	double c_2 = V.get_y() - U_.get_y();

	for (int i = 0; i < 2; i++)
	{
		double tempSlope = (i == 0) ? (-c2 / c1) : (-c_2 / c_1);
		if (is_tangent_slope(tempSlope, bound1, bound2, dir))
		{
			distance = dist(tempSlope, u, v) + dist(tempSlope, u_, v);
			if (minSumDist > distance)
			{
				minSlope = tempSlope;
				minSumDist = distance;
			}
		}
	}

	//(iii) case
	double tempSlope = (c1 + c_1) / (c2 + c_2);
	if (is_tangent_slope(tempSlope, bound1, bound2, dir) && (c1 * tempSlope + c2) * (c_1 * tempSlope + c_2) > 0)
	{
		distance = dist(tempSlope, u, v) + dist(tempSlope, u_, v);
		if (minSumDist > distance)
		{
			minSlope = tempSlope;
			minSumDist = distance;
		}
	}
	tempSlope = (c1 - c_1) / (c2 - c_2);
	if (is_tangent_slope(tempSlope, bound1, bound2, dir) && (c1 * tempSlope + c2) * (c_1 * tempSlope + c_2) < 0)
	{
		distance = dist(tempSlope, u, v) + dist(tempSlope, u_, v);
		if (minSumDist > distance)
		{
			minSlope = tempSlope;
			minSumDist = distance;
		}
	}

	*sum = minSumDist;
	return minSlope;
}
void EVENTS::compute_min_sum(void)
{
	int u=-1, u_=-1;
	double candidate;
	double startSlope = Queue[0][0]->getSlope() , minSum = std::numeric_limits<double>::infinity();;
	for (int i = 0; i < Queue.size()-1; i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			LINE* event = Queue[i][j];

			int temp_u = event->getPath(0).back();
			int temp_u_ = event->getPath(1).back();
			
				int v = event->getV();
				if (u != -1 && u_ != -1) {
					double candidateSlope = getSlopeMinSum(startSlope, event->getSlope(),rotation[i+1], v, u, u_, &candidate);
					LINE* temp = new LINE(v, candidateSlope, event->getPath(0), event->getPath(1));

					//compute sum distance for the candidate slope
					candidate += temp->getDistanceSum();
					if (candidate < minSum) {
						minSum = candidate;
						minSumLine = temp;
					}
				}
			u = temp_u, u_ = temp_u_;
			startSlope = event->getSlope();

			/*
			//structure of shortest paths changed
			if (u != temp_u || u_ != temp_u_)
			{
				if (i != 0 && j != 0) {
					int v = event->getV();
					double candidateSlope = getSlopeMinSum(startSlope, event->getSlope(), rotation[v], v, u, u_);
					LINE* minSumLine = new LINE(v, candidateSlope, event->getPath(0), event->getPath(1));

					//compute sum distance for the candidate slope
					double candidate = minSumLine->getDistanceSum();
					if (candidate < minSum)
						minSum = candidate;
				}
				u = temp_u, u_ = temp_u_;
				startSlope = event->getSlope();
			}*/
		}
	}

	minSumLine->computeEndpointWithSlope();
	printf("stop here");
}
double EVENTS::computeMinSum(void) {
	double candidateSlope[4];
	double candidateDist[4];
	LINE* prev = Queue[0][0];

	for (int i = 0; i < Queue.size(); i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			LINE* line = Queue[i][j];
			int u_s = line->getPath(0).back();
			int u_t = line->getPath(1).back();
			int v = line->getV();

			double vx = point_list[v].get_x();
			double vy = point_list[v].get_y();

			double a = point_list[u_s].get_x()- vx;
			double b = vy - point_list[u_s].get_y();
			double c = point_list[u_t].get_x() - vx;
			double d = vy - point_list[u_t].get_y();

			if (a == 0 || c == 0)
			{
				printf("error\n");
				return -1;
			}
			candidateSlope[0] = -b / a;
			candidateSlope[1] = -d / c;

			if (b + d == 0) {
				candidateSlope[2] = std::numeric_limits<double>::infinity();
				if (a + c < 0)
					candidateSlope[2] = -candidateSlope[2];
			}
			else
				candidateSlope[2] = (a + c) / (b + d);
			if (b - d == 0)
			{
				candidateSlope[3] = std::numeric_limits<double>::infinity();
				if (a - c < 0)
					candidateSlope[3] = -candidateSlope[3];
			}
			else
				candidateSlope[3] = (a - c) / (b - d);

			double slope1 = prev->getSlope();
			double slope2 = line->getSlope();

			for (int k = 0; k < 4; k++)
				candidateDist[k] = line->Dist_from_u_to_line(0,candidateSlope[k])+line->Dist_from_u_to_line(1,candidateSlope[k]);


			prev = line;
		}
	}
}