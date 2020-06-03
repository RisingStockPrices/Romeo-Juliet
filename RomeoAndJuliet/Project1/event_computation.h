#pragma once

#include<iostream>
#include<vector>

#include "LOS.h"
#define INT_MAX 100000000
using namespace std;


class EVENTS {
	int next_line_id;
	vector<vector<LOS*>> queue;
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
			queue.push_back(vector<LOS*>());
		}
		spt[0] = _spt_s;
		spt[1] = _spt_t;
	}
	vector<int> get_shortest_path()
	{
		return shortest_path;
	}
	vector<vector<LOS*>> get_queue() {
		return queue; 
	}
	void compute_path_events();
	void compute_boundary_events();
	void compute_bend_events();
	void compute_shortest_path_to_line(int i, int j);
	void sort_boundary_events();
};

bool compare_by_angle(LOS* a, LOS* b)
{
	if (a->get_path_angle() < b->get_path_angle())
		return true;
	else
		return false;
}

bool compare_angle(LINE* a, LINE* b)
{
	return (a->getAngle() < b->getAngle());
}

/* sorts boundary events in the order they appear on the polygon 
   uses the slope information stored in LOS class
   must be called only after all the path events are set*/
void EVENTS::sort_boundary_events()
{
	
	for (int i = 0; i < queue.size()-1; i++)
	{
		LOS* path = queue[i][0];
		for (int j = 0; j < queue[i].size(); j++)
		{
			LOS* boundary = queue[i][j];
			float angle = calculate_angle_between_positive(path->get_p1(), path->get_p2(), boundary->get_p1(), boundary->get_p2());
			boundary->set_path_angle(angle);
		}
		sort(queue[i].begin(), queue[i].end(), compare_by_angle);

		/*
		//set path_angle separately
		//sort by path_angle
		sort(queue[i].begin(), queue[i].end(), compare_by_angle);
		
		float angle = calculate_angle_between_positive(shortest_path[i + 1], shortest_path[i + 2], shortest_path[i], shortest_path[i + 1]);
		if (angle < queue[i].back()->get_path_angle())
		{
			vector<LOS*> sorted;
			sorted.push_back(queue[i][0]);
			
			reverse(queue[i].begin(), queue[i].end());
	
			sorted.insert(sorted.end(), queue[i].begin(), queue[i].end()-1);

			queue[i] = sorted;
		}*/
	}

	printf("sorting boundary events is complete\n");
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
		/*
		//float angle = i == 0 ? 0 : calculate_angle_between_positive(shortest_path[i], shortest_path[i + 1], shortest_path[i], shortest_path[i - 1]);
		LOS* los = new LOS(next_line_id++, prev, cur, cur, 0, PATH);
		los->compute_other_endpoint(true);
		//los->extend_path_event();
		queue[i].push_back(los);*/
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
					double angle = calculate_angle_between_positive(cur, vertex_id, prev, cur);
					BOUNDARY* boundary = new BOUNDARY(cur, vertex_id,angle);
					Queue[i - 1].push_back(boundary);
					/*
					float angle = calculate_angle_between_positive(cur, vertex_id, prev, cur);
					LOS* los = new LOS(next_line_id++, cur, vertex_id, cur, angle, j < s_size ? BOUNDARY_S : BOUNDARY_T);// BOUNDARY);
					los->compute_other_endpoint(false);
					queue[i-1].push_back(los);*/

				}
			}
		}

		for (int j = 0; j < Queue[i - 1].size(); j++)
		{
			sort(Queue[i - 1].begin(), Queue[i - 1].end(), compare_angle);
		}
	}
}

/* Sets the foot_bool, pi_s_l, pi_t_l for the given los in the queue */
/*
void EVENTS::compute_shortest_path_to_line(int i, int j)
{
	LOS* los = queue[i][j];

	if (los->get_type() == PATH)
	{
		vector<int> s_to_l(shortest_path.begin(), shortest_path.begin() + i + 1);
		vector<int> t_to_l(shortest_path.rbegin(), shortest_path.rbegin() + shortest_path.size() - i);
		los->set_foot_bool(true);
		los->set_pi_s_l(s_to_l);
		los->set_pi_t_l(t_to_l);
		return;
	}
	else //BOUNDARY CASE (_S or _T)
	{
		//let's first think about BOUNDARY_S
		int rotation = los->get_endpoint1();
		int polygon_vertex = los->get_endpoint2();
		Point other_vertex = los->get_other_endpoint();
		int p_edge_num = los->get_polygon_edge();

		vector<int> s_to_v = spt_s->retrieve_shortest_path(rotation);
		vector<int> s_to_e1 = spt_s->retrieve_shortest_path(p_edge_num);
		vector<int> s_to_e2 = spt_s->retrieve_shortest_path((p_edge_num + v_num - 2) % (v_num - 1));

		//first get the shortest path from s to the other_vertex
		vector<int> s_to_other_endpoint;
		int max_size = s_to_e1.size() <= s_to_e2.size() ? s_to_e1.size() : s_to_e2.size();
		int apex_idx = 0;
		for (apex_idx; apex_idx < max_size; apex_idx++)
		{
			if (s_to_e1[apex_idx] != s_to_e2[apex_idx])
				break;
		}
		s_to_other_endpoint.insert(s_to_other_endpoint.end(),s_to_e1.begin(), (s_to_e1.begin() + apex_idx));
		vector<int> temp1, temp2;
		temp1.insert(temp1.end(), s_to_e1.begin() + apex_idx, s_to_e1.end());
		temp2.insert(temp2.end(), s_to_e2.begin() + apex_idx, s_to_e2.end());
		get_remaining_path(temp1, temp2, los)


		for (int i = 0; i < s_to_e1.size(); i++)
		{
			if (s_to_e1[i] == s_to_e2[i])
			{
				s_to_other_endpoint.push_back(s_to_e1[i]);
			}
			else
			{
				//a function that directly computes the shortest path from foot to apex...? or chain
				break;
			}
		}
		//now get the shortest path from s to the line!!

		//S_TO_L�� v(endpoint1) �� other_endpoint���������� ���ؾ��ҵ�
		//T_TO_L�� v�� �Ƹ� endpoint1�� ��endpoint2 �����̿���
		los->set_foot_bool(false);
	}
	
}*/

void print_vector(vector<int> vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
		printf("%d ", vec[i]);
	}
	printf("\n");
}

bool is_polygon_edge(int diag)
{
	int p1 = diagonal_list[diag].get_origin();
	int p2 = diagonal_list[diag].get_dest();

	if (abs(p1 - p2) == 1)
		return true;
	if (p1 * p2 == 0 && p1 + p2 == v_num - 1)
		return true;
	else
		return false;
}
bool is_polygon_edge(int p1, int p2)
{
	if (abs(p1 - p2) == 1)
		return true;
	if (p1 * p2 == 0 && p1 + p2 == v_num - 1)
		return true;
	else
		return false;
}
int find_diagonal(int p1, int p2)
{
	for (int i = 0; i < diagonal_list.size(); i++)
	{
		if (diagonal_list[i].check_same_edge(p1, p2))
			return i;
	}

	return -1;
}

int opposite_tri(int current_tri, int diag)
{
	int* tri_cand = diagonal_list[diag].get_triangle();
	int new_tri;

	if (tri_cand[0] == current_tri)
		new_tri = tri_cand[1];
	else
		new_tri = tri_cand[0];

	if (new_tri == -1)
	{
		return -1;
	}
	else
		return new_tri;
}

Point compute_bend_event_endpoint(int p1, int p2, int rotation_vertex)
{
	Point foot = foot_of_perpendicular(rotation_vertex, point_list[p1], point_list[p2]);
	int foot_tri = point_state.find_triangle(foot);
	int vertex[2];

	int tri = choose_triangle(p2, p1, vertex);
	if (tri == foot_tri)
		return foot;	

	while (!is_polygon_edge(vertex[0],vertex[1]))
	{
		//set the new diag (vertex[0], vertex[1])
		int* diag_list = t_list[tri].get_d_list();
		int diag = -1;
		for (int i = 0; i < 3; i++)
		{
			int d = diag_list[i];
			if (d != -1 && diagonal_list[d].check_same_edge(vertex[0], vertex[1]))
			{
				diag = d;
				break;
			}
		}
		//find opposite triangle to diag
		tri = opposite_tri(tri, diag);
		if (tri == foot_tri)
			return foot;

		//getting the other endpoint
		Triangle t = t_list[tri];
		int* p_list = t.get_p_list();
		int other_p;
		for (int i = 0; i < 3; i++)
		{
			if (p_list[i] != vertex[0] && p_list[i] != vertex[1])
			{
				other_p = p_list[i];
				break;
			}
		}

		//setting vertex[0] and vertex[1] -> the next diag
		if (check_penetration(p1, p2, p2, vertex[0], other_p))
		{
			vertex[1] = other_p;
		}
		else if (check_penetration(p1, p2, p2, vertex[1], other_p))
		{
			vertex[0] = other_p;
		}
	}

	//diag should be the polygon edge
	return *get_line_intersection(p1, p2,vertex[0],vertex[1]);

}
LOS* add_bend_event(LOS* path_event, int rotation_vertex, bool first)
{
	//get the foot of perpendicular from the rotation vertex to line(p1,p2)
	int p1 = path_event->get_p1();
	int p2 = path_event->get_p2();
	Point foot = foot_of_perpendicular(rotation_vertex, point_list[p1], point_list[p2]);
	
	//if the foot is in the polygon boundary
	int valid = point_state.find_triangle(foot);
	if (valid != -1)
	{
		LOS los(-1, rotation_vertex, -1, rotation_vertex, -1, Bend);
		los.set_endpoint(0, foot);
		return &los;
		//los connects points rotation_vertex and foot
	}
	else //not inside the polygon ->we find the intersection with the polygon boundary
	{
		//guess what!? we already computed it!! it's part of the path event
		foot = path_event->get_endpoint(first);
		LOS los(-1, rotation_vertex, -1, rotation_vertex, -1, Bend);
		los.set_endpoint(0, foot);
		return &los;
	}

}

void EVENTS::compute_bend_events()
{
	//compute shortest path to line for all (path & boundary) events
	for (int i = 0; i < Queue.size(); i++)
	{
		//path events
		LINE* line = Queue[i][0];
		PATH* p = (PATH*)line;
		Point* endP = line->getEndpoints();
		pair<vector<int>, Point> res = shortest_path_line(endP[1], point_list[p->getV2()], spt[0]);
		line->setPathS(res.first);
		res = shortest_path_line(endP[0], point_list[p->getV()], spt[1]);
		line->setPathT(res.first);

		//boundary events
		for (int j = 1; j < Queue[i].size(); j++)
		{
			line = Queue[i][j];
			endP = line->getEndpoints();
			pair<vector<int>, Point> res = shortest_path_line(endP[0], endP[1], spt[0]);
			line->setPathS(res.first);
			res = shortest_path_line(endP[0], endP[1], spt[1]);
			line->setPathT(res.first);
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

			if (path.size() != prev->getPath(0).size())
				printf("stop here");
			it = find(path.begin(), path.end(), prev_u);
			//Type 2 scenario - u is gone from the sp
			if (it == path.end())
			{
				//new bend event with same rotation vertex as PREV and orthogonal to line u-u''
				BEND* bend = new BEND(prev->getV(), prev_u, prev___u);
				if (j = 0 && i > 0)
					Queue[i - 1].push_back(bend);
				else if (j == 0)
					printf("what the");
				else
					Queue[i].insert(Queue[i].begin() + j, bend);
			}
			else //u is still there
			{
				//Type 1 scenario (ii)
				if ((it + 1) != path.end() && prev__u == path.back())
				{	
					//new bend event with same rotation vertex as PREV and orthogonal to line u-u'
					BEND* bend = new BEND(prev->getV(), prev_u, prev__u);
					if (j = 0 && i > 0)
						Queue[i - 1].push_back(bend);
					else if (j == 0)
						printf("what the");
					else
						Queue[i].insert(Queue[i].begin() + j, bend);
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

	/*
	for (int i = 0; i < queue.size(); i++)
	{
		for (int j = 0; j < queue[i].size(); j++)
		{
			queue[i][j]->compute_shortest_path_to_los(shortest_path, spt);
		}
	}

	printf("done computing all the shortest paths\n");

	
	vector<int> prev = queue[0][0]->get_pi_s_l();// int prev_size = queue[0][0]->get_pi_s_l();
	for (int i = 0; i < queue.size(); i++)
	{
		for (int j = 0; j < queue[i].size(); j++)
		{
			vector<int> cur = queue[i][j]->get_pi_s_l();
			if (prev.size() != cur.size())
			{
				vector<int>* bigger = (prev.size() > cur.size()) ? &prev : &cur;
				int a = bigger->at(bigger->size() - 2);
				int b = bigger->at(bigger->size() - 1);
				if(a!=shortest_path[i] && b!=shortest_path[i])
					Point test = compute_bend_event_endpoint(a, b, shortest_path[i]);
			}
			prev = cur;
		}
	}
	//for every consecutive event... we have to see whether 
	//there is a change in the combinatorial structure of the path


	printf("done computing bend events\n");*/
	
}
