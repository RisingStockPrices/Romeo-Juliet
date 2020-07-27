#pragma once
#include "LOS.h"
#define INT_MAX 100000000
using namespace std;

LINE* minSumLine;
LINE* minMaxLine;

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
	void compute_min_max(void);
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
					BOUNDARY* boundary = new BOUNDARY(cur, vertex_id,spt,(j<s_size));
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
		bool type1 = false;
		std::sort(Queue[i].begin(), Queue[i].end(), compare_slope);
		vector<LINE*> newQueue;
		int curIdx = 0;
		for (; curIdx < Queue[i].size(); curIdx++)
		{
			if (curS == Queue[i][curIdx]->getSlope())
				break;
		}
		//this is for the tBEND_BOUNDARY_PATH events
		if (curIdx != Queue[i].size() - 1 && Queue[i][curIdx + 1]->getSlope() == curS)
			type1 = true;

		int nxtIdx = Queue[i].size();
		
		ROT rot = rotation[i+1];
		vector<LINE*>::iterator it = Queue[i].begin() + curIdx + 1;//points to place right after path event

		if (rot == CW) {
			newQueue.insert(newQueue.begin(), it, Queue[i].end());//path event 직후부터 끝까지 넣음
			newQueue.insert(newQueue.end(), Queue[i].begin(), it);//처음부터 path event 까지 넣음
			reverse(newQueue.begin(), newQueue.end());
			if (type1)
			{
				newQueue.insert(newQueue.begin() + 1, newQueue.back());
				newQueue.pop_back();
			}
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
			int idx = (j == 0) ? i : i + 1;
			int rot = shortest_path[idx];
			//inspect paths with root _s
			int idx_prev = 0, idx_cur = 0;
			while (idx_cur != cur.size() || idx_prev != prev.size())
			{
				if (idx_cur != cur.size() && idx_prev != prev.size() && cur[idx_cur] == prev[idx_prev])
					idx_cur++, idx_prev++;
				else
				{
					if (idx_cur < cur.size()) {
						vector<int>::iterator it = find(prev.begin(), prev.end(), cur[idx_cur]);
						if (it == prev.end()) //type 1 (ii)
						{
							/*
							BEND* bend = new BEND(rot, cur[idx_cur - 1], cur[idx_cur], slope_prev, slope_cur, rotation[idx]);
							Point* endP = bend->getEndpoints();
							pair<vector<int>, Point> res;
							for (int k = 0; k < 2; k++)
							{
								res = shortest_path_line(endP[0], endP[1], spt[k]);
								bend->setPath(k, res);
							}

							idx_cur++;
							
							if (j == 0)
								Queue[i - 1].push_back(bend);
							else
								tempQueue.push_back(bend);
							
							*/
							BEND* bend = new BEND(rot, cur[idx_cur - 1], cur[idx_cur], spt,0,true);
							idx_cur++;
							if (bend->getType() != tERROR && is_tangent_slope(bend->getSlope(),slope_prev,slope_cur,rotation[idx]))
							{
								if (j == 0)
									Queue[i - 1].push_back(bend);
								else
									tempQueue.push_back(bend);
							}
						}
					}
					if (idx_prev < prev.size()) {
						vector<int>::iterator it = find(cur.begin(), cur.end(), prev[idx_prev]);
						if (it == cur.end()) //type 2
						{
							/*
							BEND* bend = new BEND(rot, prev[idx_prev - 1], prev[idx_prev], slope_prev,slope_cur, rotation[idx]);
							Point* endP = bend->getEndpoints();
							pair<vector<int>, Point> res;
							for (int k = 0; k < 2; k++)
							{
								res = shortest_path_line(endP[0], endP[1], spt[k]);
								bend->setPath(k, res);
							}

							idx_prev++;
							if (j == 0)
								Queue[i - 1].push_back(bend);
							else
								tempQueue.push_back(bend);
							*/
							BEND* bend = new BEND(rot, prev[idx_prev - 1], prev[idx_prev], spt,0,false);
							idx_prev++;
							if (bend->getType() != tERROR && is_tangent_slope(bend->getSlope(),slope_prev,slope_cur,rotation[idx]))
							{
								if (j == 0)
									Queue[i - 1].push_back(bend);
								else
									tempQueue.push_back(bend);
							}
						}
					}

				}
			}
			//bend events in terms of paths with root _t
			idx_prev = 0, idx_cur = 0;
			while (idx_cur != cur_t.size() || idx_prev != prev_t.size())
			{
				if (idx_cur != cur_t.size() && idx_prev != prev_t.size() && cur_t[idx_cur] == prev_t[idx_prev])
					idx_cur++,idx_prev++;
				else
				{
					if (idx_cur < cur_t.size()) {
						vector<int>::iterator it = find(prev_t.begin(), prev_t.end(), cur_t[idx_cur]);
						if (it == prev_t.end()) //type 1 (ii)
						{
							BEND* bend = new BEND(rot, cur_t[idx_cur - 1], cur_t[idx_cur], spt,1,false);
							idx_cur++;
							if (bend->getType() != tERROR && is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[idx]))
							{
								if (j == 0)
									Queue[i - 1].push_back(bend);
								else
									tempQueue.push_back(bend);
							}
						}
					}
					if (idx_prev < prev_t.size()) {
						vector<int>::iterator it = find(cur_t.begin(), cur_t.end(), prev_t[idx_prev]);
						if (it == cur_t.end()) //type 2
						{
							BEND* bend = new BEND(rot, prev_t[idx_prev - 1], prev_t[idx_prev], spt,1,true);
							idx_prev++;
							if (bend->getType() != tERROR && is_tangent_slope(bend->getSlope(), slope_prev, slope_cur, rotation[idx]))
							{
								if (j == 0)
									Queue[i - 1].push_back(bend);
								else
									tempQueue.push_back(bend);
							}
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

	
	ROT before = DEFAULT;
	for (int i = 1; i < Queue.size(); i++)
	{
		ROT current = rotation[i];
		if (before!=DEFAULT && before != current)
		{
			Queue[i - 1][0]->setType1(true);
			/*
			PATH* path = (PATH*) Queue[i-1][0];
			BEND* type1 = new BEND(path);
			vector<int> route = shortest_path_random_point(point_list[path->getV()], spt[0]);
			type1->setPath(0, route);
			route = shortest_path_random_point(point_list[path->getV2()], spt[1]);
			type1->setPath(1, route);
			Queue[i - 1].push_back(type1);*/
		}
		before = current;
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


double getSlopeMinSum(double bound1, double bound2, ROT dir, int v, int u, int u_, double* sum)
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
double getSlopeMinMax(LINE* prev, LINE* cur, ROT dir, int v, int u, int u_, double* Max)
{
	double bound1 = prev->getSlope();
	double bound2 = prev->getSlope();

	double minDist = max(prev->getLength(0), prev->getLength(1));
	double minSlope;

	Point U = point_list[u];
	Point U_ = point_list[u_];
	Point V = point_list[v];
	double c1 = U.get_x() - V.get_x();
	double c2 = V.get_y() - U.get_y();
	double c_1 = U_.get_x() - V.get_x();
	double c_2 = V.get_y() - U_.get_y();

	//check common point9
	double Bsquare = (prev->getLength_noFoot(0) - prev->getLength_noFoot(1));
	Bsquare *= Bsquare;
	double D1 = (c1 + c_1), D2 = (c2 + c_2);
	double a = D1*D1-Bsquare;
	double b = 2*D1*D2;
	double c = D2 * D2 - Bsquare;

	double D = b * b - 4 * a*c;
	double tempSlope;
	
	if (D == 0)
		tempSlope = -b / (2*a);
	else if(D>0){
		double first = (-b + sqrt(D)) / (2*a);
		double second = (-b - sqrt(D)) / (2*a);

		if ((c1*first + c2)*(c_1*first + c_2) < 0)
			tempSlope = first;
		else if ((c1*second + c2)*(c_1*second + c_2) < 0)
			tempSlope = second;
		else {
			printf("this shouldn't be happening\n");
			exit(34);
		}
	}

	if (D>=0 && is_tangent_slope(tempSlope, bound1, bound2, dir))
	{
		*Max = min(dist(tempSlope, u, v) + prev->getLength(0), dist(tempSlope, u_, v)+prev->getLength(1));
		return tempSlope;
	}
	else
	{
		//가장 자리 확인 띠
		double first = max(dist(bound1, u, v) + prev->getLength(0), dist(bound1, u_, v) + prev->getLength(1));
		double second = max(dist(bound2, u, v) + cur->getLength(0), dist(bound2, u_, v) + cur->getLength(1));
		
		*Max = min(first, second);
		return (first > second) ? second : first;
	}
}
double getSlopeMinSum(LINE * prev, LINE * cur, ROT dir, int v, int u, int u_, double* sum)
{
	double bound1 = prev->getSlope();
	double bound2 = cur->getSlope();
	//case (i) : boundaries
	double minDist = min(prev->getLength(), cur->getLength());
	double minSlope;// = (minDist == bound1) ? bound1 : bound2;
	
	Point U = point_list[u];
	Point U_ = point_list[u_];
	Point V = point_list[v];
	double c1 = U.get_x() - V.get_x();
	double c2 = V.get_y() - U.get_y();
	double c_1 = U_.get_x() - V.get_x();
	double c_2 = V.get_y() - U_.get_y();

	//case (ii) : absolute value boundaries
	for (int i = 0; i < 2; i++)
	{
		double tempSlope = (i == 0) ? (-c2 / c1) : (-c_2 / c_1);
		if (is_tangent_slope(tempSlope, bound1, bound2, dir))
		{
			double distance = dist(tempSlope, u, v) + dist(tempSlope, u_, v);
			if (minDist > distance)
				minSlope = tempSlope, minDist = distance;
		}
	}

	//case (iii) local optimum/minimum candidates
	double tempSlope = (c1 + c_1) / (c2 + c_2);
	if (is_tangent_slope(tempSlope, bound1, bound2, dir) && (c1 * tempSlope + c2) * (c_1 * tempSlope + c_2) > 0)
	{
		double distance = dist(tempSlope, u, v) + dist(tempSlope, u_, v);
		if (minDist > distance)
		{
			minSlope = tempSlope;
			minDist = distance;
		}
	}
	tempSlope = (c1 - c_1) / (c2 - c_2);
	if (is_tangent_slope(tempSlope, bound1, bound2, dir) && (c1 * tempSlope + c2) * (c_1 * tempSlope + c_2) < 0)
	{
		double distance = dist(tempSlope, u, v) + dist(tempSlope, u_, v);
		if (minDist > distance)
		{
			minSlope = tempSlope;
			minDist = distance;
		}
	}

	*sum = minDist;
	return minSlope;
}
void EVENTS::compute_min_max(void) {
	//find v
	int v, i;
	double slope_prev, slope_next;
	for (i = 0; i < Queue.size()-1; i++)
	{
		LINE* line = Queue[i][0];
		if (line->getLength(0) >= line->getLength(1)) {
			v = line->getV();
			slope_prev = line->getSlope();
			slope_next = Queue[i + 1][0]->getSlope();
			break;
		}
	}

	LINE* start = Queue[i][0], *end = Queue[i][0], *prev = Queue[i][0];
	double minMax = std::numeric_limits<double>::infinity();
	vector<double> tempSum;
	LINE* reference;
	for (int j = 0; j < Queue[i].size(); j++)
	{
		end = Queue[i][j];
		TYPE type_start = start->getType();
		TYPE type_end = end->getType();

		if (type_end == tBOUNDARY)
			continue;

		int idx = (j == 0) ? i : i + 1;
		if (type_start == tBEND_del || start->getType1())
			reference = end;
		else
			reference = start;

		int u = reference->getPath(0).back();
		int u_ = reference->getPath(1).back();
		int v = reference->getV();
		if (type_end == tPATH)
			v = start->getV();

		double sum;

		//needs to be changed to min max
		double slope = getSlopeMinMax(start, end, rotation[idx], v, u, u_, &sum);
		//double slope = getSlopeMinSum(start, cur, rotation[idx], refV, refU, refU_, &sum);
		LINE * temp = new LINE(v, slope, reference->getPath(0), reference->getPath(1)); //need to adjust this later for bend_del cases
		sum = max(temp->getLength(0)+dist(slope,u,v), temp->getLength(1)+dist(slope,u_,v));// temp->getDistanceSum();
		//tempSum.push_back(sum);
		if (sum < minMax) {
			minMax = sum; minMaxLine = temp;
		}

		start = end;
	}

}
void EVENTS::compute_min_sum(void)
{
	LINE* start = Queue[0][0], * end = Queue[0][0], * prev = Queue[0][0];
	double minSum = std::numeric_limits<double>::infinity();
	vector<double> tempSum;
	LINE* reference;
	for (int i = 0; i < Queue.size(); i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			end = Queue[i][j];
			TYPE type_start = start->getType();
			TYPE type_end = end->getType();

			if (type_end == tBOUNDARY)
				continue;

			int idx = (j == 0) ? i : i + 1;
			if (type_start == tBEND_del || start->getType1())
				reference = end;
			else
				reference = start;
			
			int u = reference->getPath(0).back();
			int u_ = reference->getPath(1).back();
			int v = reference->getV();
			if (type_end == tPATH)
				v = start->getV();

			double sum;
			double slope = getSlopeMinSum(start, end, rotation[idx], v, u, u_, &sum);
			//double slope = getSlopeMinSum(start, cur, rotation[idx], refV, refU, refU_, &sum);
			LINE * temp = new LINE(v, slope, reference->getPath(0), reference->getPath(1)); //need to adjust this later for bend_del cases
			sum += temp->getDistanceSum();
			tempSum.push_back(sum);
			if (sum < minSum) {
				minSum = sum; minSumLine = temp;
			}

			start = end;

			/*
			cur = Queue[i][j];
			idx = (j == 0) ? i : i + 1;
			v = cur->getV();
			u = cur->getPath(0).back();
			u_ = cur->getPath(1).back();

			prev_u = (prev->getType() == tBEND_del) ? u:prev_u;
			prev_u_ = (prev->getType() == tBEND_del) ? u_ : prev_u_;
			bool exception = (cur->getType() == tBEND_del);
			
			bool change = prev_u != u || prev_u_ != u_ || prev_v != v;
			if(exception || change)
			{
				bool back = (start->getType() == tBEND_del);
				LINE* reference;
				if (back)
				{
					reference = cur;
				}
				else
					reference = start;

				int refV = reference->getV();
				int refU = reference->getPath(0).back();
				int refU_ = reference->getPath(1).back();

				double sum;
				double slope = getSlopeMinSum(start, cur, rotation[idx], refV, refU, refU_, &sum);
				LINE * temp = new LINE(reference->getV(), slope, reference->getPath(0), reference->getPath(1)); //need to adjust this later for bend_del cases
				sum += temp->getDistanceSum();
				if (sum < minSum) {
					minSum = sum; minSumLine = temp;
				}

				start = cur;
				prev_v = v;
				prev_u = u;
				prev_u_ = u_;
			}

			prev = cur;
		}
	}



	
	for (int i = 0; i < Queue.size(); i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			next = Queue[i][j];
			idx = (j == 0) ? i : i + 1;

			/*
			// compute u, u_, v
			v = next->getV();
			u = next->getPath(0).back();
			u_ = next->getPath(1).back();
			

			if (v != prev_v || u != prev_u || u_ != prev_u_)
			{
				int final_v = prev_v;
				int final_u = prev_u;
				int final_u_ = prev_u_;

				//do the chunk computation for start~end
				start;
				end;

				//set start and end for next chunk
				start = cur;
				prev_v = v, prev_u = u, prev_u_ = u_;
			}
			
			end = next;
			cur = next;
			prev = cur;
			*/

			/*
			v = cur->getV();
			u = cur->getPath(0).back();
			u_ = cur->getPath(1).back();
			if (cur->getType() == tBEND_del || start->getType1()) {
				u = next->getPath(0).back(), u_ = next->getPath(1).back();
			}

			if (prev_u != u || prev_u_ != u_ || prev_v != v)
			{
				//compute for start ~end
				double sum;
				double slope = getSlopeMinSum(start, end, rotation[idx], prev_v, prev_u, prev_u_, &sum);
				LINE* reference = (start->getType() == tBEND_del || start->getType1()) ? end : start;
				LINE * temp = new LINE(prev_v, slope, reference->getPath(0), reference->getPath(1)); //need to adjust this later for bend_del cases
				sum += temp->getDistanceSum();
				if (sum < minSum) {
					minSum = sum; minSumLine = temp;
				}

				start = cur, end = next;
				prev_u = u, prev_u_ = u_, prev_v = v;
			}
			else
			{
				end = next;
			}
			cur = next;
			prev = cur;*/
		}
	}
	/*
	//don't forget the final case [ CUR should be the last path event ]
	v = cur->getV(), u = cur->getPath(0).back(), u_ = cur->getPath(1).back();
	if (prev_u != u || prev_u_ != u_ || prev_v != v)
	{
		double sum;
		double slope = getSlopeMinSum(start, end, rotation[Queue.size() - 1], prev_v, prev_u, prev_u_, &sum);
		LINE * reference = (start->getType() == tBEND_del || start->getType1()) ? end : start;
		LINE * temp = new LINE(prev_v, slope, reference->getPath(0), reference->getPath(1));
		sum += temp->getDistanceSum();
		if (sum < minSum)
			minSum = sum, minSumLine = temp;
	}*/

	//LINE* prev = Queue[0][0], *cur = Queue[0][0];
	//LINE* literal_prev = prev;
	//int cur_u, cur_u_, cur_v;
	//int prev_u = prev->getPath(0).back(), prev_u_ = prev->getPath(1).back(), prev_v = prev->getV();
	/*
	bool was_bend = false, isType1 = false;
	double candidate;
	double startSlope = Queue[0][0]->getSlope() , minSum = std::numeric_limits<double>::infinity();

	LINE* next,*start=Queue[0][0],*cur = Queue[0][0];

	int start_u = start->getPath(0).back(), start_v = start->getV(), start_u_ = start->getPath(1).back();

	for (int i = 0; i < Queue.size(); i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			int next_i = (j == Queue[i].size() - 1) ? i + 1 : i;
			int next_j = (j == Queue[i].size() - 1) ? 0 : j + 1;
			int idx = (j == 0) ? i : i + 1;
			next = Queue[next_i][next_j];
			int next_u = next->getPath(0).back();
			int next_u_ = next->getPath(1).back();
			int next_v = next->getV();

			int v = (next_j == 0) ? start_v : next_v;
			if (next_u != start_u || next_u_ != start_u_ || v != start_v)
			{
				int u, u_, v;
				double candidateSlope = getSlopeMinSum(start, cur, rotation[idx], start_u, start_u_, start_v, &candidate);
				LINE* temp = new LINE(shortest_path[idx], candidateSlope,cur->getPath(0), cur->getPath(1));
				candidate += temp->getDistanceSum();
				if (candidate < minSum)
					minSum = candidate, minSumLine = temp;
				start = cur;
				start_u = next_u;
				start_v = (cur->getType()==tPATH)?start_v:next_v;
				start_u_ = next_u_;
			}

			cur = next;
		}
	}
	for (int i = 0; i < Queue.size(); i++)
	{
		for (int j = 0; j < Queue[i].size(); j++)
		{
			cur = Queue[i][j];
			cur_u = cur->getPath(0).back();
			cur_u_ = cur->getPath(1).back();
			cur_v = cur->getV();
			int idx = (j == 0) ? i : i + 1;
			if (cur_v == prev_v && cur_u == prev_u && cur_u_ == prev_u_) {
				literal_prev = cur;
				continue;
			}
			else {
				int u = prev_u;
				int u_ = prev_u_;
				int v = prev_v;//shortest_path[idx];
				double candidateSlope = getSlopeMinSum(prev,literal_prev, rotation[idx], u, u_, v,&candidate);
				LINE* temp = new LINE(shortest_path[idx], candidateSlope, cur->getPath(0), cur->getPath(1));
				candidate += temp->getDistanceSum();
				if (candidate < minSum)
					minSum = candidate, minSumLine = temp;
				prev = cur;
				prev_u = cur_u, prev_u_ = cur_u_, prev_v = cur_v;
				literal_prev = cur;
			}
		}
	}*/


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

			double a = point_list[u_s].get_x() - vx;
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
				candidateDist[k] = line->Dist_from_u_to_line(0, candidateSlope[k]) + line->Dist_from_u_to_line(1, candidateSlope[k]);


			prev = line;
		}
	}
}