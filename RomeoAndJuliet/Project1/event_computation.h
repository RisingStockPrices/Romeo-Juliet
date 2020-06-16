#pragma once

#include<iostream>
#include<vector>
#include <algorithm>
#include "LOS.h"
#define INT_MAX 100000000
using namespace std;

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
		
		PATH* line = new PATH(prev, cur);
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
					BOUNDARY* boundary = new BOUNDARY(cur, vertex_id);
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
	bool inBetween = (direction == CW) == (from > to);
	if (inBetween)
		return (from - slope) * (to - slope) < 0;
	else
		return (from - slope) * (to - slope) > 0;
}
void EVENTS::compute_bend_events()
{
	//compute shortest path to line for all (path & boundary) events
	for (int i = 0; i < Queue.size(); i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			LINE* line = Queue[i][j];
			Point* endP = line->getEndpoints();
			pair<vector<int>, Point> res = shortest_path_line(endP[0], endP[1], spt[0]);
			line->setPath(0, res);
			res = shortest_path_line(endP[0], endP[1], spt[1]);
			line->setPath(1, res);
		}
	}

	LINE* prev;
	int prev_u = -1;
	int prev__u = -1;
	int prevu = Queue[0][0]->getPath(0).back();
	
	for (int i = 0; i < Queue.size(); i++)
	{
		vector<LINE*> BendEvents;
		for (int j = 0; j < Queue[i].size(); j++)
		{
			LINE* cur = Queue[i][j];
			vector<int> curPath = cur->getPath(0);

			//last vertex on path to line
			int u = curPath.back();
			//vertex following u in shortest path
			vector<int>::iterator it = find(shortest_path.begin(), shortest_path.end(), u);
			int _u = (it != shortest_path.end() && (it+1)!=shortest_path.end()) ? *(it + 1) : -1;
			int __u = (curPath.size() > 2) ? curPath[curPath.size() - 2]:-1;
			
			//possible bend events
			if (prevu != u)
			{
				//type 2 if u==prev__u
				//type 1 (ii) if u==prev_u
				int rot = prev->getV();
				int prevIdx = (j == 0) ? i - 1 : i;
				
				//making sure not type 1 (i) case (since it's already added)
				if (rotation[prevIdx] == rotation[prevIdx + 1])
				{
					BEND* bend = new BEND(rot, u, prevu, 0);
					if (bend->getType() != tERROR)
					{
						//check for tangent
						bool isTangent = is_tangent_slope(bend->getSlope(), prev->getSlope(), cur->getSlope(), rotation[prevIdx + 1]);
						if (isTangent)
						{
							if (j == 0)
								Queue[prevIdx].push_back(bend);
							else
								BendEvents.push_back(bend);//Queue[prevIdx].push_back(bend);
						}
					}
				}
			}
			prevu = u;
			prev_u = _u;
			prev__u = __u;
			prev = cur;
		}
		Queue[i].insert(Queue[i].end(),BendEvents.begin(),BendEvents.end());
	}

	prevu = Queue.back().back()->getPath(1).back();
	int s_size = shortest_path.size();
	for (int i = Queue.size() - 1; i >= 0; i--)
	{
		vector<LINE*> BendEvents;
		for (int j = Queue[i].size() - 1; j >= 0; j--)
		{
			LINE* cur = Queue[i][j];
			if (cur->getType() > tBOUNDARY)
				continue;

			vector<int> curPath = cur->getPath(1);
			//last vertex on path to line
			int u = curPath.back();
			//vertex following u in shortest path
			vector<int>::iterator it = find(shortest_path.begin(), shortest_path.end(), u);
			int _u = (it != shortest_path.end() && it!=shortest_path.begin()) ? *(it - 1) : -1;
			int __u = (curPath.size() > 2) ? curPath[curPath.size() - 2] : -1;

			//possible bend events
			if (prevu != u)
			{
				//type 2 if u==prev__u
				//type 1 (ii) if u==prev_u
				int rot = prev->getV();

				//making sure not type 1 (i) case (since it's already added)
 				if (i+3<s_size && rotation[i+1] == rotation[i+ 2])
				{
					BEND* bend = new BEND(rot, u, prevu, 0);
					if (bend->getType() != tERROR)
					{
						//check for tangent
						bool isTangent = is_tangent_slope(bend->getSlope(), prev->getSlope(), cur->getSlope(), (ROT)(3-rotation[i + 1]));
						if (isTangent)
							BendEvents.push_back(bend);
					}
				}
			}

			prevu = u;
			prev_u = _u;
			prev__u = __u;
			prev = cur;
		}
		Queue[i].insert(Queue[i].end(), BendEvents.begin(), BendEvents.end());
	}

	/*
	//u and _u each correspond to u and u' in the paper (see page 7 - bend events)
	int u, _u;
	LINE* prev = Queue[0][0];
	vector<int> prev_path = prev->getPath(0);
	int prev_u = prev_path.back();
	vector<int>::iterator it = find(shortest_path.begin(), shortest_path.end(), prev_u);
	int prev__u = *(it + 1);
	int prev___u = -1;

	//inspecting shortest paths from S
	for (int i = 0; i < Queue.size(); i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			LINE* line = Queue[i][j];
			vector<int> path = line->getPath(0);
			bool isBendEvent = false;

			int rotIdx = distance(shortest_path.begin(), it);
			it = find(path.begin(), path.end(), prev_u);
			int orth2 = -1;
			if (it == path.end())//type 2 (spline vertex deleted)
				orth2 = prev___u;
			else if ((it + 1) != path.end() && prev__u == path.back()) //type 1 (ii) : a new vertex added in the spline - the vertex should follow u in the sp
			{
				orth2 = prev__u;
				
				//check rotation vertex's rotation direction see if it mismatches previous in path
			}
			
			if (orth2!=-1)
			{
				int prev_v = prev->getV();
				BEND* bend = new BEND(prev_v, prev_u, orth2,0);
				if (bend->getType() != tERROR) {
					Point endP = bend->getEndpoints()[0];
					point_list.push_back(endP);
					if (i + 2 >= shortest_path.size())
						printf("start by admitting from cradle to tomb\n");
					else {
						if (j == 0 && i > 0)
						{
							bool isTangent = is_tangent(shortest_path[i - 1], shortest_path[i], shortest_path[i + 1], point_list.size() - 1);
							if (isTangent) {
								Queue[i - 1].push_back(bend);
							}
						}
						else
						{
							bool isTangent = is_tangent(shortest_path[i], shortest_path[i + 1], shortest_path[i + 2], point_list.size() - 1);
							if (isTangent) {
								Queue[i].push_back(bend);// insert(Queue[i].begin() + j, bend);
							}
						}
					}
					point_list.pop_back();
				}
				else
				{
					printf("wait\n");
				}
			}

			prev_u = (path.empty()) ? -1 : path.back();
			prev___u = (path.size() > 1) ? *(path.end() - 2) : -1;
			it = find(shortest_path.begin(), shortest_path.end(), prev_u);
			if (it == shortest_path.end() || it + 1 == shortest_path.end())
				prev__u = -1;
			else
				prev__u = *(it + 1);
			prev = line;
		}
	}

	prev = Queue[Queue.size()-1][0];
	prev_path = prev->getPath(1);
	prev_u = prev_path.back();
	it = find(shortest_path.begin(), shortest_path.end(), prev_u);
	int sp_size = shortest_path.size();
	prev__u = *(it - 1);
	prev___u = -1;

	for (int i = Queue.size()-1;i>=0;i--)//0; i < Queue.size(); i++)
	{
		for (int j = Queue[i].size()-1; j >= 0; j--)//; j < Queue[i].size(); j++)
		{
			LINE* line = Queue[i][j];
			if (line->getType() != tBEND) {
				vector<int> path = line->getPath(1);
				bool isBendEvent = false;

				it = find(path.begin(), path.end(), prev_u);
				int orth2 = -1;
				if (it == path.end())//type 2 (spline vertex deleted)
				{
					orth2 = prev___u;
					isBendEvent = true;
				}
				else if ((it + 1) != path.end() && prev__u == path.back())
				{
					orth2 = prev__u;
					isBendEvent = true;
				}

				if (isBendEvent)
				{
					int prev_v = prev->getV();
					BEND* bend = new BEND(prev_v, prev_u, orth2, 0);
					if (bend->getType() != tERROR) {
						Point endP = bend->getEndpoints()[0];
						point_list.push_back(endP);
						if (i + 2 >= shortest_path.size())
							printf("start by admitting from cradle to tomb\n");
						else {

							bool isTangent = is_tangent(shortest_path[i], shortest_path[i+1], shortest_path[i + 2], point_list.size() - 1);
							
							if (isTangent)
								Queue[i].push_back(bend);// insert(Queue[i].begin() + j, bend);						
						}
						point_list.pop_back();
					}
				}

				prev_u = (path.empty()) ? -1 : path.back();
				prev___u = (path.size() > 1) ? *(path.end() - 2) : -1;
				it = find(shortest_path.begin(), shortest_path.end(), prev_u);
				if (it == shortest_path.end() || it == shortest_path.begin())
					prev__u = -1;
				else
					prev__u = *(it - 1);
				prev = line;
			}
		}
	}
	*/
	sort_by_slope();
}

double EVENTS::computeMinSum(void)
{
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