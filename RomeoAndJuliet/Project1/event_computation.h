#pragma once

#include<iostream>
#include<vector>

#include "LOS.h"
#define INT_MAX 100000000
using namespace std;


class EVENTS {
	int next_line_id;
	vector<vector<LINE*>> Queue;
	vector<int> shortest_path;
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
	void compute_path_events();
	void compute_boundary_events();
	void compute_bend_events();
};


bool compare_angle(LINE* a, LINE* b)
{
	return (a->getAngle() < b->getAngle());
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
					double angle = 0;
					if (j >= s_size)
						angle = calculate_angle_between(cur, vertex_id, cur, prev);

					else
						angle = calculate_angle_between(cur, vertex_id, prev, cur);
					angle = abs(normalize_angle(angle));
					BOUNDARY* boundary = new BOUNDARY(cur, vertex_id,angle);
					Queue[i - 1].push_back(boundary);
				}
			}
		}

		for (int j = 0; j < Queue[i - 1].size(); j++)
			sort(Queue[i - 1].begin(), Queue[i - 1].end(), compare_angle);
		
	}
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

	//u and _u each correspond to u and u' in the paper (see page 7 - bend events)
	int u, _u;
	LINE* prev = Queue[0][0];
	vector<int> prev_path = prev->getPath(0);
	int prev_u = prev_path.back();
	vector<int>::iterator it = find(shortest_path.begin(), shortest_path.end(), prev_u);
	int prev__u = *(it + 1);
	int prev___u = -1;

	for (int i = 0; i < Queue.size(); i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			LINE* line = Queue[i][j];
			vector<int> path = line->getPath(0);
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
							if (isTangent)
								Queue[i - 1].push_back(bend);
						}
						else
						{
							bool isTangent = is_tangent(shortest_path[i], shortest_path[i + 1], shortest_path[i + 2], point_list.size() - 1);
							if(isTangent)
								Queue[i].insert(Queue[i].begin() + j, bend);

						}
					}
					point_list.pop_back();
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
							{
				
								if (j == 0 && i > 0)
									Queue[i - 1].push_back(bend);
								else		
									Queue[i].insert(Queue[i].begin() + j, bend);
							}
						
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
	
}
