#pragma once
#include "Hourglass.h"
#include <algorithm>
#include <set>
#include "Point.h"
#include "Edge.h"
#include "polygon_decomposition.h"
#include "DebugPrint.h"
#include "tangent.h"

SNode * find_middle_diagonal(SNode * root, int e1, int e2) {
	if (root->get_diagonal() == -1)
		return root;

	SNode ** children = root->get_children();
	bool check = true;
	for (int i = 0; i < 2; i++) {
		check = children[i]->check_inclusive(e1);
		check = check && children[i]->check_inclusive(e2);
		if (check) {
			return find_middle_diagonal(children[i], e1, e2);
		}
	}
	return root;
}
String * duplicate_string(String *s, Chain ** c, int apax) {
	vector<int> p_1 = c[0]->get_point_list();
	vector<int> p_2 = c[1]->get_point_list();

	int it1, it2, diff_it1, diff_it2;
	if (apax == p_1.front()) {
		it1 = 0;
		diff_it1 = 1;
	}
	else {
		it1 = p_1.size() - 1;
		diff_it1 = -1;
	}

	if (apax == p_2.front()) {
		it2 = 0;
		diff_it2 = 1;
	}
	else {
		it2 = p_2.size() - 1;
		diff_it2 = -1;
	}

	it1 = it1 + diff_it1;
	it2 = it2 + diff_it2;
	while (it1 >= 0 && it1 < (int)p_1.size() && it2 >= 0 && it2 < (int)p_2.size()) {
		if (p_1[it1] != p_2[it2]) {
			break;
		}
		it1 = it1 + diff_it1;
		it2 = it2 + diff_it2;
	}
	it1 = it1 - diff_it1;
	it2 = it2 - diff_it2;
	apax = p_1[it1];
	delete(c[0]);
	vector<int> new_s;
	if (diff_it1 == 1) {
		c[0] = new Chain(vector<int>(p_1.begin() + it1, p_1.end()));
		new_s = vector<int>(p_1.begin(), p_1.begin() + it1 + 1);
	}
	else {
		int t1 = p_1.size() - it1;
		c[0] = new Chain(vector<int>(p_1.begin(), p_1.end() - t1 + 1));
		new_s = vector<int>(p_1.end() - t1, p_1.end());
	}
	delete(c[1]);
	if (diff_it2 == 1) {
		c[1] = new Chain(vector<int>(p_2.begin() + it2, p_2.end()));
	}
	else {
		int t2 = p_2.size() - it2;
		c[1] = new Chain(vector<int>(p_2.begin(), p_2.end() - t2 + 1));
	}
	String * ret = NULL;
	if (new_s.size() >  1) {
		if (s == NULL)
			ret = new String(new_s);
		else ret = ret = new String(s, new String(new_s));
	}
	return ret;
}
void Hourglass::duplicate_strings() {
	String * new_s = NULL;
	if (first_chain[0] != NULL) {
		if (s == NULL) {
			int f00 = first_chain[0]->get_point(0), f01 = first_chain[0]->get_last_point();
			int f10 = first_chain[1]->get_point(0), f11 = first_chain[1]->get_last_point();
			int f_a, s_a = -1;
			if (edge_list[0].check_same_point(f00) == -1 && (f00 == f10 || f00 == f11)) {
				s_a = f00;
				f_a = f00;
				new_s = duplicate_string(s, first_chain, f_a);
			}
			else if (edge_list[0].check_same_point(f01) == -1 && (f01 == f10 || f01 == f11)) {
				f_a = f01;
				s_a = f01;
				new_s = duplicate_string(s, first_chain, f_a);
			}

			if (new_s != NULL) {
				s = new_s;
				apex[0] = f_a;
				apex[1] = s_a;
				second_chain[0] = new Chain(s_a);
				second_chain[1] = new Chain(s_a);
			}
		}
		else {
			new_s = duplicate_string(s, first_chain, apex[0]);
			if (new_s != NULL) {
				s = new_s;
			}
		}
	}
	if (second_chain[0] != NULL) new_s = duplicate_string(s, second_chain, apex[1]);
	if (new_s != NULL) s = new_s;

	return;
}
vector<Hourglass> hourglass_list;
Hourglass concatenate_hourglasses(Hourglass& left, Hourglass& right);
Hourglass concatenate_hourglasses(int h1, int h2);
Hourglass construct_hourglass_point(int p, int high) {
	Edge e(p);
	diagonal_list.push_back(e);
	//hourglass_list.push_back(Hourglass(diagonal_list.size()-1, high));
	return Hourglass(diagonal_list.size() - 1, high);//��� �߰��� point edge �� funnel? ����� ��
}
int construct_hourglass(int low, int high) {
	low %= v_num;
	high %= v_num;
	SNode * node = diagonal_list[low].get_SNode();
	//-1 return -> triangle
	SNode * middle = find_middle_diagonal(node, low, high);
	Hourglass new_h;
	int middle_diagonal = middle->get_diagonal();
	if (middle_diagonal == -1) {
		new_h = Hourglass(low, high, middle);
	}
	else {
		int h1 = s_graph[middle_diagonal][low];
		if (h1 == -2) h1 = construct_hourglass(middle_diagonal, low);

		int h2 = s_graph[middle_diagonal][high];
		if (h2 == -2) h2 = construct_hourglass(middle_diagonal, high);

		new_h = concatenate_hourglasses(h1, h2);
	}
	new_h.set_id();
	hourglass_list.push_back(new_h);
	s_graph[low][high] = new_h.get_id();
	return new_h.get_id();
}
void init_hourglass_val() {
	id_num = 0;
	hourglass_list = vector<Hourglass>();
}
String* concatenate_two_funnels_cc(Chain** left, int apax1, Chain** right, int apax2, Edge * common_edge) {

	Hourglass h_left;
	h_left.set_first_edge(Edge(apax1));
	h_left.set_second_edge(*common_edge);
	h_left.set_first_chain(left);

	Hourglass h_right;
	h_right.set_first_edge(Edge(apax2));
	h_right.set_second_edge(*common_edge);
	h_right.set_first_chain(right);

	Edge * temp = common_edge;
	Hourglass new_h = concatenate_hourglasses(h_right, h_left);
	return new_h.get_string();
	/*String * min = NULL;
	int t1, t2;
	for(int i=0;i<2;i++)
	for (int j = 0; j < 2; j++) {
	compute_tangent(left[i], right[j],&t1,&t2);
	vector<int> v_left = left[i]->get_point_list();
	vector<int> v_right = right[j]->get_point_list();
	vector<int> new_v;
	if (v_left[0] != common_edge->get_origin() && v_left[0] != common_edge->get_dest()) {
	new_v.insert(new_v.end(), v_left.begin(), v_left.begin() + t1 +1);
	}
	else {
	new_v.insert(new_v.end(), v_left.begin()+t1,v_left.end());
	reverse(new_v.begin(), new_v.end());
	}

	if (v_right[0] != common_edge->get_origin() && v_right[0] != common_edge->get_dest()) {
	if (new_v.back() == *(v_right.begin() + t2)) {
	vector<int> temp(v_right.begin(), v_right.begin() + t2);
	reverse(temp.begin(), temp.end());
	new_v.insert(new_v.end(), temp.begin(), temp.end());
	}
	else {
	vector<int> temp(v_right.begin(), v_right.begin() + t2+1);
	reverse(temp.begin(), temp.end());
	new_v.insert(new_v.end(), temp.begin(), temp.end());
	}
	}
	else {
	if (new_v.back() == *(v_right.begin() + t2)) {
	new_v.insert(new_v.end(), v_right.begin() + t2, v_right.end());
	}
	else {
	new_v.insert(new_v.end(), v_right.begin() + t2 + 1, v_right.end());
	}
	}

	String * new_string = new String(new_v);

	if (min == NULL || min->get_length() > new_string->get_length())
	min = new_string;
	else delete(new_string);
	}
	return min;*/
}
Return_val concatenate_two_funnels_oc(Chain** left, Chain** right, Edge * common_edge, int apax, Edge right_edge) {
	Hourglass h_left;
	h_left.set_first_edge(Edge(apax));
	h_left.set_second_edge(*common_edge);
	h_left.set_first_chain(left);

	Hourglass h_right;
	h_right.set_first_edge(right_edge);
	h_right.set_second_edge(*common_edge);
	h_right.set_first_chain(right);

	Edge * temp = common_edge;
	Hourglass new_h = concatenate_hourglasses(h_right, h_left);
	common_edge = temp;
	return Return_val(new_h.get_string(), new_h.get_apaxes()[0], new_h.get_first_chain());
}
Return_val concatenate_two_funnels_oc(Chain** left, Chain** right, Edge * common_edge) {
	String * min = NULL;
	int t1, t2;
	String * s_list[2][2];
	Chain * c_list[2][2][2];
	int apax[2][2];
	bool valid[2][2] = { true, true, true, true };
	for (int i = 0; i<2; i++)
		for (int j = 0; j < 2; j++) {
			vector<int> v_left = left[i]->get_point_list();
			vector<int> v_right = right[j]->get_point_list();
			int v_l0 = v_left.front(), v_l1 = v_left.back(), v_r0 = v_right.front(), v_r1 = v_right.back();
			//outer
			if (v_l0 == v_r0 || v_l0 == v_r1 || v_l1 == v_r0 || v_l1 == v_r1) {
				bool check = compute_outer_tangent(left[i], right[j], &t1, &t2, common_edge, left[(i + 1) % 2], right[(j + 1) % 2]);
				if (check == false) {
					valid[i][j] = false;
					continue;
				}
			}//inner
			else {
				bool check = compute_inner_tangent(left[i], right[j], &t1, &t2, common_edge, left[(i + 1) % 2], right[(j + 1) % 2]);
				if (check == false) {
					valid[i][j] = false;
					continue;
				}
			}

			vector<int> new_v;
			if (common_edge->check_same_point(v_left[0]) == -1) {
				new_v.insert(new_v.end(), v_left.begin(), v_left.begin() + t1 + 1);
			}
			else {
				new_v.insert(new_v.end(), v_left.begin() + t1, v_left.end());
				reverse(new_v.begin(), new_v.end());
			}
			s_list[i][j] = new String(new_v);

			apax[i][j] = v_left[t1];
			c_list[i][j][0] = right[j]->cutting_chain(common_edge, apax[i][j], right[j]);
			c_list[i][j][1] = right[(j + 1) % 2]->cutting_chain(common_edge, apax[i][j], right[(j + 1) % 2]);
		}
	int min_i, min_j;
	point_type min_val = -1;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			if (valid[i][j] == false) continue;
			String * s = duplicate_string(s_list[i][j], c_list[i][j], apax[i][j]);
			if (s != NULL) s_list[i][j] = s;
			point_type len_val = 2 * s_list[i][j]->get_length() + c_list[i][j][0]->get_len() + c_list[i][j][1]->get_len();
			if (min_val == -1 || min_val > len_val) {
				min_val = len_val;
				min_i = i;
				min_j = j;
			}
		}
	}
	return Return_val(s_list[min_i][min_j], apax[min_i][min_j], c_list[min_i][min_j]);
}
Hourglass concatenate_hourglasses(int h1, int h2) {
	return concatenate_hourglasses(hourglass_list[h1], hourglass_list[h2]);
}



Edge computeOuterTangent(vector<int> left_point_list, vector<int> right_point_list, vector<int> sum, int common, bool upper) {
	int samePoint = common;
	bool point_tangent = true;

	bool(*one_sided_test)(vector<int>, int, int) = upper ? all_left : all_right;

	if (sum.size() == 2)
		return Edge(sum[0], sum[1]);

	for (int i = 0; i < sum.size() - 1; i++)
	{
		if (!one_sided_test(sum, sum[i], sum[i + 1]))
		{
			point_tangent = false;
			break;
		}
	}

	if (point_tangent)
	{
		return Edge(samePoint, samePoint);
	}
	else {
		for (int i = 0; i < left_point_list.size(); i++) {
			for (int j = 0; j < right_point_list.size(); j++) {
				if (one_sided_test(sum, left_point_list[i], right_point_list[j])) {
					if (left_point_list[i] == right_point_list[j])
						samePoint = right_point_list[j];
					else
						return Edge(left_point_list[i], right_point_list[j]);
				}
			}
		}

		return Edge(samePoint, samePoint);
	}
}

enum state {
	neutral,
	pending,
	upward,
	downward,
};

bool set_current_state(state& current_state, point_type x, point_type y, Point p1, Point p2)
{
	current_state = neutral;

	point_type p1_y = p1.get_y();
	point_type p2_y = p2.get_y();
	point_type p1_x = p1.get_x();
	point_type p2_x = p2.get_x();

	if (p1.get_y() == y && p1_x > x)
	{
		if (p1_y > p2_y)
		{
			current_state = downward;
		}
		else if (p1_y < p2_y)
		{
			current_state = upward;
		}
		else
			current_state = pending;

		return true;
	}
	if (p2.get_y() == y && p2_x>x)
	{
		if (p1_y > p2_y)
		{
			current_state = upward;
		}
		else if (p1_y < p2_y)
		{
			current_state = downward;
		}
		else
			current_state = pending;
		return true;
	}
	else {
		return false;
	}

}

bool check_inclusive(vector<int> chain, int test_point)
{
	Point test = point_list[test_point];
	Point horizontal(10000, test.get_y());//MAX_x + 1, test.get_y()); //do forgive me...

	if (chain.size() == 2)//���ϼ��� test_point�� ���� ���� true �� return �ؾߴ�!!
	{
		if (check_line_intersection_open(test, test, point_list[chain[0]], point_list[chain[1]]))
			return true;
		else
			return false;
	}
	chain.push_back(chain.front());

	point_type test_x = test.get_x();
	point_type test_y = test.get_y();

	int count = 0;

	Point from = point_list[chain[0]];
	Point to;

	state previous_state = neutral;
	state current_state = neutral;

	for (int i = 0; i < chain.size() - 1; i++)
	{
		to = point_list[chain[i + 1]];

		if (set_current_state(current_state, test_x, test_y, from, to))
		{
			if (previous_state == upward)
			{
				previous_state = neutral;
				if (current_state == pending)
					previous_state = upward;
				else if (current_state == downward)
					count++;
			}
			else if (previous_state == downward)
			{
				previous_state = neutral;
				if (current_state == pending)
					previous_state = downward;
				else if (current_state == upward)
					count++;
			}
			else
			{
				if (previous_state == neutral && current_state == pending)
					count++;
				previous_state = current_state;
			}
		}
		else {//current state�� neutral �̱� ������ previous_state�� ���� ������!
			if (check_line_intersection_open(test, horizontal, from, to))//check_line_intersection(test, horizontal, from, to))
				count++;
		}
		from = to;
	}
	if (count % 2 == 0)//false
		return false;
	else//true
		return true;
}
int check_inclusive(vector<int> chain, int test_point, int neighbor)
{
	Point test = point_list[test_point];
	Point horizontal(10000, test.get_y());//MAX_x + 1, test.get_y()); //do forgive me...
	int edge_num = -1;

	if (chain.size() == 2)//���ϼ��� test_point�� ���� ���� true �� return �ؾߴ�!!
	{
		if (check_line_intersection_open(test, test, point_list[chain[0]], point_list[chain[1]]))
			return 0;
		else
			return -1;
	}

	chain.push_back(chain.front());

	point_type test_x = test.get_x();
	point_type test_y = test.get_y();

	int count = 0;

	Point from = point_list[chain[0]];
	Point to;

	state previous_state = neutral;
	state current_state = neutral;

	for (int i = 0; i < chain.size() - 1; i++)
	{
		to = point_list[chain[i + 1]];

		if (neighbor == -1)
			edge_num = 0;
		else if (i<chain.size() - 2 && check_line_intersection(test_point, neighbor, chain[i], chain[i + 1],true))
			edge_num = i;
		if (set_current_state(current_state, test_x, test_y, from, to))
		{
			if (previous_state == upward)
			{
				previous_state = neutral;
				if (current_state == pending)
					previous_state = upward;
				else if (current_state == downward)
					count++;
			}
			else if (previous_state == downward)
			{
				previous_state = neutral;
				if (current_state == pending)
					previous_state = downward;
				else if (current_state == upward)
					count++;
			}
			else
			{
				if (previous_state == neutral && current_state == pending)
					count++;
				previous_state = current_state;
			}
		}
		else {//current state�� neutral �̱� ������ previous_state�� ���� ������!
			if (check_line_intersection_open(test,horizontal,from,to))//check_line_intersection(test, horizontal, from, to))
			{
				count++;
			}
		}
		from = to;
	}
	if (count % 2 == 0)//false
		return -1;
	else//true
	{
		if (edge_num == -1)
			exit(11);
		return edge_num;
	}
}

vector<int> cross_the_line_check(Edge tangent, vector<int> other_side, vector<int> same_side, int common_point, bool upper)
{
	bool(*which_side)(int, int, int) = upper ? is_left : is_right;
	vector<int> hump;
	vector<int> pointList;

	if (tangent.is_point())
	{
		hump = vector<int>(); //is this necessary?
		return hump;
	}
	int from = (tangent.is_reverse()) ? tangent.get_dest() : tangent.get_origin();
	int to = (tangent.is_reverse()) ? tangent.get_origin() : tangent.get_dest();
	//find sub vector starting & ending with tangent
	vector<int>::iterator from_it = find(same_side.begin(), same_side.end(), from);
	vector<int>::iterator to_it = find(same_side.begin(), same_side.end(), to);
	//they should enter this function already ordered (left->right)
	//set up pointList (bounding area)
	pointList.insert(pointList.end(), from_it, to_it);
	pointList.push_back(to);

	for (int i = 0; i < other_side.size(); i++)
	{
		if (find(pointList.begin(), pointList.end(), other_side[i]) == pointList.end()) {
			if (check_inclusive(pointList, other_side[i]))//don't call when testpoint is in chain
				hump.push_back(other_side[i]);
		}
	}

	return hump;
}

vector<int> mountanizeHump(Edge tangent, vector<int> outliers, vector<int> same_side, bool upper)
{
	vector<int> temp;
	vector<int> boundary;
	vector<int> mountain;
	vector<pair<int, int>> same_side_outliers;

	int edge_num;
	int start = tangent.is_reverse() ? tangent.get_dest() : tangent.get_origin();
	int end = tangent.is_reverse() ? tangent.get_origin() : tangent.get_dest();

	bool(*valid_tangent_test)(vector<int>, int, int) = upper ? all_right : all_left;
	temp.push_back(start);
	temp.insert(temp.end(), outliers.begin(), outliers.end());
	temp.push_back(end);

	vector<int>::iterator from = temp.begin();
	vector<int>::iterator to;

	//set the boundary first, considering only the tangent & outliers from the other side (not the same_side)
	boundary.push_back(start);
	while (*from != end)
	{
		for (to = from + 1; to != temp.end(); to++)
		{
			if (valid_tangent_test(temp, *from, *to))//!one_side_test(temp, *from, *to, upper))
			{
				boundary.push_back(*to);
				from = to;
				break;
			}
		}
	}//boundary set

	int start_idx = -1, end_idx = -1;
	//need to get the outliers from the same side that invade the boundary
	for (int i = 0; i < same_side.size(); i++)
	{
		if (same_side[i] == start)
			start_idx = i;
		else if (same_side[i] == end)
			end_idx = i;
	}

	for (int idx = start_idx + 1; idx < end_idx; idx++) //tangent ����(?) �� �ִ� ������ Ȯ���ϸ� �Ǵµ�
	{
		edge_num = check_inclusive(boundary, same_side[idx], same_side[idx - 1]);

		if (edge_num != -1)
			same_side_outliers.push_back(make_pair(edge_num, same_side[idx]));
	}

	int boundary_index = -1;
	for (int i = 0; i < same_side_outliers.size(); i++)
	{
		edge_num = same_side_outliers[i].first;
		for (int j = boundary_index + 1; j <= edge_num; j++)
			mountain.push_back(boundary[j]);

		mountain.push_back(same_side_outliers[i].second);
		boundary_index = edge_num;
	}

	//������ �κе鵵 push!!
	for (int i = boundary_index + 1; i < boundary.size(); i++)
	{
		mountain.push_back(boundary[i]);
	}
	return mountain;
}
Chain* invalid_outer_chains(Edge tangent, vector<int> left_chain, vector<int> right_chain, vector<int> same_side, vector<int> outliers, Edge leftEdge, Edge rightEdge, bool upper)
{
	int left_tangent_point = tangent.is_reverse() ? tangent.get_dest() : tangent.get_origin();
	int right_tangent_point = tangent.is_reverse() ? tangent.get_origin() : tangent.get_dest();

	vector<int> piAB;
	vector<int> tanBC;
	vector<int> piCD;

	Chain* result = new Chain();

	vector<int>::iterator it = find(left_chain.begin(), left_chain.end(), left_tangent_point);
	if (it == left_chain.end())
	{
		exit(1);
	}

	if (left_chain[0] == leftEdge.get_origin() || left_chain[0] == leftEdge.get_dest())
	{
		piAB.insert(piAB.begin(), left_chain.begin(), it);
	}
	else {
		piAB.insert(piAB.begin(), it, left_chain.end());
		reverse(piAB.begin(), piAB.end());
		piAB.pop_back();
		//reverse(piAB.begin(), piAB.end());
	}

	it = find(right_chain.begin(), right_chain.end(), right_tangent_point);
	if (it == right_chain.end())
	{
		exit(2);
	}

	if (right_chain[0] == rightEdge.get_origin() || right_chain[0] == rightEdge.get_dest())
	{
		piCD.insert(piCD.begin(), right_chain.begin(), it);
		reverse(piCD.begin(), piCD.end());
	}
	else
	{
		piCD.insert(piCD.begin(), it + 1, right_chain.end());
	}

	tanBC = mountanizeHump(tangent, outliers, same_side, upper);

	result->append_points(piAB);
	result->append_points(tanBC);
	result->append_points(piCD);
	return result;
}
Chain* valid_outer_chains(Edge tangent, vector<int> left_chain, vector<int> right_chain, Edge leftEdge, Edge rightEdge)
{
	int left_tangent_point = tangent.is_reverse() ? tangent.get_dest() : tangent.get_origin();
	int right_tangent_point = tangent.is_reverse() ? tangent.get_origin() : tangent.get_dest();

	vector<int> piAB;
	vector<int> tanBC;
	vector<int> piCD;

	Chain* result = new Chain();

	vector<int>::iterator it = find(left_chain.begin(), left_chain.end(), left_tangent_point);
	if (it == left_chain.end())
	{
		exit(1);
	}

	if (left_chain[0] == leftEdge.get_origin() || left_chain[0] == leftEdge.get_dest())
	{
		piAB.insert(piAB.begin(), left_chain.begin(), it);
	}
	else {
		piAB.insert(piAB.begin(), it, left_chain.end());
		reverse(piAB.begin(), piAB.end());
		piAB.pop_back();
	}

	it = find(right_chain.begin(), right_chain.end(), right_tangent_point);
	if (it == right_chain.end())
	{
		exit(2);
	}

	if (right_chain[0] == rightEdge.get_origin() || right_chain[0] == rightEdge.get_dest())
	{
		piCD.insert(piCD.begin(), right_chain.begin(), it);
		reverse(piCD.begin(), piCD.end());
	}
	else
	{
		piCD.insert(piCD.begin(), it + 1, right_chain.end());
	}

	tanBC.push_back(left_tangent_point);
	if (!tangent.is_point())
		tanBC.push_back(right_tangent_point);

	result->append_points(piAB);
	result->append_points(tanBC);
	result->append_points(piCD);
	return result;
}

int nearest_point_to_common(Chain* chain, Edge common)
{
	vector<int> points = chain->get_point_list();

	//check origin
	if (points.front() == common.get_origin() || points.front() == common.get_dest())
		return points[1];//second to first
	else
		return points[points.size() - 2];//second to last
}

bool check_sharks_fin_case(Edge common, Edge other)
{
	if (is_right(other.get_origin(), common.get_origin(), common.get_dest()) && is_right(other.get_dest(), common.get_origin(), common.get_dest()))
	{
		return false;
	}
	else if (is_left(other.get_origin(), common.get_origin(), common.get_dest()) && is_left(other.get_dest(), common.get_origin(), common.get_dest()))
	{
		return false;
	}
	else
	{
		return true;
	}
}
Hourglass concatenateOpenOpen(Hourglass& _left, Hourglass& _right)
{
	//set first_edge and second edge
	Hourglass newHourglass;

	Edge* leftEdgeList = _left.get_edge_list();
	Edge* rightEdgeList = _right.get_edge_list();
	Chain** leftChainList = _left.get_first_chain();//first_chain���� ���� ������ �Ŵϱ� �̷��� �ӽô�.
	Chain** rightChainList = _right.get_first_chain();

	Edge commonEdge;
	int leftCommonEdgeIndex = 0, rightCommonEdgeIndex = 0;

	for (int i = 0; i < 2; i++) //�� edge ��� ���� ���� ���ٰ� �����Դϴ�.
	{
		for (int j = 0; j < 2; j++)
		{
			if (leftEdgeList[i] == rightEdgeList[j]) {
				leftCommonEdgeIndex = i;
				rightCommonEdgeIndex = j;
				commonEdge = Edge(leftEdgeList[i].get_origin(), leftEdgeList[i].get_dest());
				break;
			}
		}
	}

	Edge leftEdge = leftEdgeList[!leftCommonEdgeIndex];
	Edge rightEdge = rightEdgeList[!rightCommonEdgeIndex];

	if (leftEdge == rightEdge)
		return _left;

	//left hourglass�� common���� ���� edge�� first edge�� set��!! -> �� ������ �߿��Ѱ�?

	newHourglass.set_first_edge(leftEdge);
	newHourglass.set_second_edge(rightEdge);

	bool left_sharks_fin = check_sharks_fin_case(commonEdge, leftEdge);
	bool right_sharks_fin = check_sharks_fin_case(commonEdge, rightEdge);

	int common_upper_point = commonEdge.get_origin();
	int common_lower_point = commonEdge.get_dest(); //if this doesn't work, we'll swap the two

	int which_chain_to_test = 0;
	bool found = false;
	bool left = true;
	Chain* chain_to_check;

	Chain* left_upper_chain;
	Chain* left_lower_chain;
	Chain* right_upper_chain;
	Chain* right_lower_chain;

	if (!left_sharks_fin && !right_sharks_fin)
	{
		if (leftEdge.is_point() && commonEdge.check_same_point(leftEdge.get_dest()) != -1)//the special case �̰� true �̸� right�� �׷� �� ����!
		{
			which_chain_to_test = 2;
		}
		while (!found)
		{
			switch (which_chain_to_test) {
			case 0:
			case 1:
				chain_to_check = leftChainList[which_chain_to_test];
				break;
			case 2:
			case 3:
				left = false;
				chain_to_check = rightChainList[which_chain_to_test - 2];
				break;
			default:
				printf("this shouldn't be happening\n");
				exit(10);
			}

			if (chain_to_check->get_point_list_size() != 1)
				found = true;
			else
				which_chain_to_test++;
		}

		int nearest = nearest_point_to_common(chain_to_check, commonEdge);
		bool(*side)(int, int, int) = left ? is_left : is_right;

		if (!side(nearest, common_lower_point, common_upper_point))
			swap(common_upper_point, common_lower_point);

		//the common Edge's upper and lower point should be correctly set.....

	}
	else//�� �� �ϳ��� ������ �϶�
	{
		if (left_sharks_fin)
		{
			Chain* chain_with_origin = leftChainList[0]->check_inclusive(leftEdge.get_origin()) ? leftChainList[0] : leftChainList[1];
			if (is_right(leftEdge.get_dest(), leftEdge.get_origin(), commonEdge.get_origin()) && is_right(leftEdge.get_dest(), leftEdge.get_origin(), commonEdge.get_dest()))
			{
				if (chain_with_origin->check_inclusive(commonEdge.get_origin()))
				{
					common_upper_point = commonEdge.get_origin();
					common_lower_point = commonEdge.get_dest();
				}
				else
				{
					common_upper_point = commonEdge.get_dest();
					common_lower_point = commonEdge.get_origin();
				}
			}
			else
			{
				if (chain_with_origin->check_inclusive(commonEdge.get_origin()))
				{
					common_upper_point = commonEdge.get_dest();
					common_lower_point = commonEdge.get_origin();
				}
				else
				{
					common_upper_point = commonEdge.get_origin();
					common_lower_point = commonEdge.get_dest();
				}
			}
		}
		else
		{
			Chain* chain_with_origin = rightChainList[0]->check_inclusive(rightEdge.get_origin()) ? rightChainList[0] : rightChainList[1];
			if (is_left(rightEdge.get_dest(), rightEdge.get_origin(), commonEdge.get_origin()) && is_left(rightEdge.get_dest(), rightEdge.get_origin(), commonEdge.get_dest()))
			{
				if (chain_with_origin->check_inclusive(commonEdge.get_origin()))
				{
					common_upper_point = commonEdge.get_origin();
					common_lower_point = commonEdge.get_dest();
				}
				else
				{
					common_upper_point = commonEdge.get_dest();
					common_lower_point = commonEdge.get_origin();
				}
			}
			else
			{
				if (chain_with_origin->check_inclusive(commonEdge.get_origin()))
				{
					common_upper_point = commonEdge.get_dest();
					common_lower_point = commonEdge.get_origin();
				}
				else
				{
					common_upper_point = commonEdge.get_origin();
					common_lower_point = commonEdge.get_dest();
				}
			}
		}
	}


	if (leftChainList[0]->check_inclusive(common_upper_point) && leftChainList[1]->check_inclusive(common_lower_point))
	{
		left_upper_chain = leftChainList[0];
		left_lower_chain = leftChainList[1];
	}
	else
	{
		left_upper_chain = leftChainList[1];
		left_lower_chain = leftChainList[0];
	}
	if (rightChainList[0]->check_inclusive(common_upper_point) && rightChainList[1]->check_inclusive(common_lower_point))
	{
		right_upper_chain = rightChainList[0];
		right_lower_chain = rightChainList[1];
	}
	else
	{
		right_upper_chain = rightChainList[1];
		right_lower_chain = rightChainList[0];
	}

	/*
	int left_chain_with_common_origin = leftChainList[0]->check_inclusive(commonEdge.get_origin()) ? 0 : 1;
	int right_chain_with_common_origin = rightChainList[0]->check_inclusive(commonEdge.get_origin()) ? 0 : 1;
	int left_chain_with_left_origin = leftChainList[0]->check_inclusive(leftEdge.get_origin()) ? 0 : 1;
	int right_chain_with_right_origin = rightChainList[0]->check_inclusive(rightEdge.get_origin()) ? 0 : 1;
	int common_upper_point = -1;
	int common_lower_point = -1;;
	int upper_index, lower_index;
	if (leftEdge.is_point())
	{
	lower_index = is_left(commonEdge.get_dest(), leftEdge.get_dest(), commonEdge.get_origin()) ? left_chain_with_common_origin : !left_chain_with_common_origin;
	}
	else
	{
	lower_index = is_left(leftEdge.get_dest(), leftEdge.get_origin(), commonEdge.get_origin()) && is_left(leftEdge.get_dest(), leftEdge.get_origin(), commonEdge.get_dest()) ? left_chain_with_left_origin : !left_chain_with_left_origin;
	}
	upper_index = !lower_index;

	left_upper_chain = leftChainList[upper_index];
	left_lower_chain = leftChainList[lower_index];

	if (rightEdge.is_point())
	{
	upper_index = is_left(commonEdge.get_dest(), rightEdge.get_dest(), commonEdge.get_origin()) ? right_chain_with_common_origin : !right_chain_with_common_origin;
	}
	else
	{
	upper_index = is_left(rightEdge.get_dest(), rightEdge.get_origin(), commonEdge.get_origin()) && is_left(rightEdge.get_dest(), rightEdge.get_origin(), commonEdge.get_dest()) ? right_chain_with_right_origin : !right_chain_with_right_origin;
	}
	lower_index = !upper_index;

	right_upper_chain = rightChainList[upper_index];
	right_lower_chain = rightChainList[lower_index];

	bool leftEdge_common_point=leftEdge.is_point() && commonEdge.check_same_point(leftEdge.get_origin())!=-1;
	bool rightEdge_common_point= rightEdge.is_point() && commonEdge.check_same_point(rightEdge.get_origin())!=-1;
	if (!leftEdge_common_point)
	{
	common_upper_point = left_upper_chain->check_inclusive(commonEdge.get_origin()) ? commonEdge.get_origin() : commonEdge.get_dest();
	common_lower_point = commonEdge.get_dest() == common_upper_point ? commonEdge.get_origin() : commonEdge.get_dest();
	}
	else if (!rightEdge_common_point)
	{
	common_upper_point = right_upper_chain->check_inclusive(commonEdge.get_origin()) ? commonEdge.get_origin() : commonEdge.get_dest();
	common_lower_point = commonEdge.get_dest() == common_upper_point ? commonEdge.get_origin() : commonEdge.get_dest();
	}
	else
	{
	exit(8);
	//not sure what this case is gonna look like
	}

	if (leftEdge_common_point)
	{
	if (leftChainList[0]->check_inclusive(common_upper_point) && leftChainList[1]->check_inclusive(common_lower_point))
	{
	upper_index = 0;
	}
	else
	{
	upper_index = 1;
	}

	lower_index = !upper_index;
	left_upper_chain = leftChainList[upper_index];
	left_lower_chain = leftChainList[lower_index];
	}
	if (rightEdge_common_point)
	{
	if (rightChainList[0]->check_inclusive(common_upper_point) && rightChainList[1]->check_inclusive(common_lower_point))
	{
	upper_index = 0;
	}
	else
	{
	upper_index = 1;
	}

	lower_index = !upper_index;
	right_upper_chain = rightChainList[upper_index];
	right_lower_chain = rightChainList[lower_index];
	}*/

	/*
	//choose common_upper and common_lower
	if (!leftEdge.is_point())
	{
	int lower_index = is_left(leftEdge.get_dest(), leftEdge.get_origin(), commonEdge.get_origin()) && is_left(leftEdge.get_dest(), leftEdge.get_origin(), commonEdge.get_dest()) ? left_chain_with_left_origin : !left_chain_with_left_origin;
	int upper_index = !lower_index;

	left_lower_chain = leftChainList[lower_index];
	left_upper_chain = leftChainList[upper_index];

	common_upper_point = left_upper_chain->check_inclusive(commonEdge.get_origin()) ? commonEdge.get_origin() : commonEdge.get_dest();
	common_lower_point = commonEdge.get_dest() == common_upper_point ? commonEdge.get_origin() : commonEdge.get_dest();

	if (rightEdge.is_point())
	{

	}
	else
	{

	}
	}
	else//leftEdge �� point�� rightEdge�� ���ؾ��Ѵ�
	{
	if (!rightEdge.is_point())
	{
	int upper_index = is_left(rightEdge.get_dest(), rightEdge.get_origin(), commonEdge.get_origin()) && is_left(rightEdge.get_dest(), rightEdge.get_origin(), commonEdge.get_dest()) ? right_chain_with_right_origin : !right_chain_with_right_origin;
	int lower_index = !upper_index;

	right_upper_chain = rightChainList[upper_index];
	right_lower_chain = rightChainList[lower_index];

	common_upper_point = right_upper_chain->check_inclusive(commonEdge.get_origin()) ? commonEdge.get_origin() : commonEdge.get_dest();
	common_lower_point = commonEdge.get_dest() == common_upper_point ? commonEdge.get_origin() : commonEdge.get_dest();

	}
	else
	{

	}
	}*/

	/*
	if (leftEdge.is_point())
	{
	if (is_left(commonEdge.get_dest(), leftEdge.get_dest(), commonEdge.get_origin())) {//�̰� is_right���� �ʿ����??
	left_lower_chain = leftChainList[left_chain_with_common_origin];
	left_upper_chain = leftChainList[!left_chain_with_common_origin];
	}
	else {
	left_upper_chain = leftChainList[left_chain_with_common_origin];
	left_lower_chain = leftChainList[!left_chain_with_common_origin];
	}
	}
	else {
	if (is_left(leftEdge.get_dest(), leftEdge.get_origin(), commonEdge.get_origin()) && is_left(leftEdge.get_dest(), leftEdge.get_origin(), commonEdge.get_dest())) {
	left_lower_chain = leftChainList[left_chain_with_left_origin];
	left_upper_chain = leftChainList[!left_chain_with_left_origin];
	}
	else {
	left_upper_chain = leftChainList[left_chain_with_left_origin];
	left_lower_chain = leftChainList[!left_chain_with_left_origin];
	}
	}
	if (rightEdge.is_point())
	{
	if (is_left(commonEdge.get_dest(), rightEdge.get_dest(), commonEdge.get_origin())) {
	right_upper_chain = rightChainList[right_chain_with_common_origin];
	right_lower_chain = rightChainList[!right_chain_with_common_origin];
	}
	else {
	right_lower_chain = rightChainList[right_chain_with_common_origin];
	right_upper_chain = rightChainList[!right_chain_with_common_origin];
	}
	}
	else {
	if (is_left(rightEdge.get_dest(), rightEdge.get_origin(), commonEdge.get_origin()) && is_left(rightEdge.get_dest(), rightEdge.get_origin(), commonEdge.get_dest())) {
	right_upper_chain = rightChainList[right_chain_with_right_origin];
	right_lower_chain = rightChainList[!right_chain_with_right_origin];
	}
	else {
	right_lower_chain = rightChainList[right_chain_with_right_origin];
	right_upper_chain = rightChainList[!right_chain_with_right_origin];
	}
	}*/


	vector<int> left_upper_points = left_upper_chain->get_point_list();
	vector<int> left_lower_points = left_lower_chain->get_point_list();
	vector<int> right_lower_points = right_lower_chain->get_point_list();
	vector<int> right_upper_points = right_upper_chain->get_point_list();
	//assign left&right / upper&lower chains ^^^^^^^^^^

	/*
	for (int i = 0; i < left_upper_points.size(); i++)
	{
	if (left_upper_points[i] == v_num)
	left_upper_points[i] = 0;;
	}
	for (int i = 0; i < left_lower_points.size(); i++)
	{
	if (left_lower_points[i] == v_num)
	left_lower_points[i] =0;
	//left_lower_points[i] = left_lower_points[i] % v_num;
	}
	for (int i = 0; i < right_upper_points.size(); i++)
	{
	if (right_upper_points[i] == v_num)
	right_upper_points[i] = 0;
	//right_upper_points[i] = right_upper_points[i] % v_num;
	}
	for (int i = 0; i < right_lower_points.size(); i++)
	{
	if (right_lower_points[i] == v_num)
	right_lower_points[i] =0;
	//right_lower_points[i] = right_lower_points[i] % v_num;
	}*/


	//connect them!
	vector<int> upper_points;
	vector<int> lower_points;
	for (int i = 0; i < upper_points.size(); i++)
	{
		printf("%d ", upper_points[i]);
	}
	cout << endl;
	for (int i = 0; i < lower_points.size(); i++)
	{
		printf("%d ", lower_points[i]);
	}
	cout << endl;

	int check = connect_vectors(left_upper_points, right_upper_points, upper_points);
	if (check == -1)
	{
		printf("what the...");
	}
	check = connect_vectors(left_lower_points, right_lower_points, lower_points);
	if (check == -1)
	{
		printf("gotcha!");
	}

	Edge upperT = computeOuterTangent(left_upper_points, right_upper_points, upper_points, common_upper_point, true);
	Edge lowerT = computeOuterTangent(left_lower_points, right_lower_points, lower_points, common_lower_point, false);


	//computed outerTangents -> not all valid yet
	vector<int> UpperOutliers = cross_the_line_check(upperT, lower_points, upper_points, common_lower_point, true);
	vector<int> LowerOutliers = cross_the_line_check(lowerT, upper_points, lower_points, common_upper_point, false);

	Chain* upperChain = new Chain();
	Chain* lowerChain = new Chain();

	Chain* first[2];
	Chain* second[2];

	bool open = true;
	if (UpperOutliers.empty())
	{
		upperChain = valid_outer_chains(upperT, left_upper_points, right_upper_points, leftEdge, rightEdge);
	}
	else {
		open = false;
		upperChain = invalid_outer_chains(upperT, left_upper_points, right_upper_points, upper_points, UpperOutliers, leftEdge, rightEdge, true);
	}

	if (LowerOutliers.empty())
	{
		lowerChain = valid_outer_chains(lowerT, left_lower_points, right_lower_points, leftEdge, rightEdge);
	}
	else {
		open = false;
		lowerChain = invalid_outer_chains(lowerT, left_lower_points, right_lower_points, lower_points, LowerOutliers, leftEdge, rightEdge, false);
	}

	if (upperChain->get_point_list() == lowerChain->get_point_list()) //for the point-point edge case -> should be closed
		open = false;

	vector<int> string;
	std::set<int> upperS;
	std::set<int> lowerS;
	vector<int> final_upper_list = upperChain->get_point_list();
	vector<int> final_lower_list = lowerChain->get_point_list();
	vector<int>::iterator up;
	vector<int>::iterator low;
	vector<int> upper;
	vector<int> lower;
	int start, end;

	if (open)//open case~
	{
		first[0] = upperChain;
		first[1] = lowerChain;
		newHourglass.set_first_chain(first);
	}
	else {
		up = final_upper_list.begin();
		low = final_lower_list.begin();
		lowerS.insert(*low);

		while (up != final_upper_list.end() || low != final_lower_list.end())
		{
			if (up != final_upper_list.end()) {
				if (lowerS.find(*up) != lowerS.end())
				{
					start = *up;
					break;
				}
				else
					upperS.insert(*up);
			}
			if (low != final_lower_list.end()) {
				if (upperS.find(*low) != upperS.end())
				{
					start = *low;
					break;
				}
				else
					lowerS.insert(*low);
			}

			if (up != final_upper_list.end())
				up++;
			if (low != final_lower_list.end())
				low++;
		}//setting the start of the string

		up = find(final_upper_list.begin(), final_upper_list.end(), start);
		low = find(final_lower_list.begin(), final_lower_list.end(), start);


		upper.insert(upper.begin(), final_upper_list.begin(), up);
		lower.insert(lower.begin(), final_lower_list.begin(), low);
		upper.push_back(start);
		lower.push_back(start);
		first[0] = new Chain(upper);
		first[1] = new Chain(lower);

		while (up != final_upper_list.end() && low != final_lower_list.end() && *up == *low)
		{
			end = *up;
			string.push_back(end);
			up++;
			low++;
		}

		upper.clear();
		lower.clear();
		upper.push_back(end);
		lower.push_back(end);
		upper.insert(upper.end(), up, final_upper_list.end());
		lower.insert(lower.end(), low, final_lower_list.end());
		second[0] = new Chain(upper);
		second[1] = new Chain(lower);

		newHourglass.set_string(new String(string));
		newHourglass.set_first_chain(first);
		newHourglass.set_second_chain(second);
		newHourglass.set_apex(start, 0);
		newHourglass.set_apex(end, 1);
	}

	printSummary(_left, _right, newHourglass, open);

	return newHourglass;
	//return newHourglass;
}

Hourglass concatenate_hourglasses(Hourglass& _left, Hourglass& _right) {

	bool left_openess = _left.check_openess();
	bool right_openess = _right.check_openess();
	Edge* left_edge_list = _left.get_edge_list();
	Edge* right_edge_list = _right.get_edge_list();

	int common_edge_check[2];

	for (int i = 0; i < 2; i++) {//�� hourglass�� �����ϴ� edge ã�Ƽ� common_edge_check�� index����
		for (int j = 0; j < 2; j++) {
			if (left_edge_list[i] == right_edge_list[j]) {
				common_edge_check[0] = i;
				common_edge_check[1] = j;
				common_edge = &left_edge_list[i];
			}
		}
	}

	//printf("%d %d\n", common_edge_check[0], common_edge_check[1]);
	Hourglass new_hourglass;

	if (!left_openess && !right_openess) {
		//set left
		Chain ** lc, **rc;
		int la, ra;
		if (common_edge_check[0] == 0) {
			new_hourglass.set_first_edge(left_edge_list[1]);
			new_hourglass.set_first_chain(_left.get_second_chain());
			new_hourglass.set_apex(_left.get_apaxes()[1], 0);
			lc = _left.get_first_chain();
			la = _left.get_apaxes()[0];
		}
		else {
			new_hourglass.set_first_edge(left_edge_list[0]);
			new_hourglass.set_first_chain(_left.get_first_chain());
			new_hourglass.set_apex(_left.get_apaxes()[0], 0);
			lc = _left.get_second_chain();
			la = _left.get_apaxes()[1];
		}
		//set right
		if (common_edge_check[1] == 0) {
			new_hourglass.set_second_edge(right_edge_list[1]);
			new_hourglass.set_second_chain(_right.get_second_chain());
			new_hourglass.set_apex(_right.get_apaxes()[1], 1);
			rc = _right.get_first_chain();
			ra = _right.get_apaxes()[0];
		}
		else {
			new_hourglass.set_second_edge(right_edge_list[0]);
			new_hourglass.set_second_chain(_right.get_first_chain());
			new_hourglass.set_apex(_right.get_apaxes()[0], 1);
			rc = _right.get_second_chain();
			ra = _right.get_apaxes()[1];
		}

		String* middle_string = concatenate_two_funnels_cc(lc, la, rc, ra, common_edge);
		String* new_string;
		if (middle_string != NULL) {
			new_string = new String(_left.get_string(), middle_string);
			new_string = new String(new_string, _right.get_string());
			new_hourglass.set_string(new_string);
		}
		else {
			new_string = new String(_left.get_string(), _right.get_string());
			new_hourglass.set_string(new_string);
		}
	}//left : close, right : open
	else if (!left_openess && right_openess) {

		Return_val r_val;
		new_hourglass.set_string(_left.get_string());
		if (common_edge_check[1] == 0) {
			new_hourglass.set_second_edge(right_edge_list[1]);
		}
		else {
			new_hourglass.set_second_edge(right_edge_list[0]);
		}
		//set left
		if (common_edge_check[0] == 0) {
			new_hourglass.set_first_edge(left_edge_list[1]);
			new_hourglass.set_first_chain(_left.get_second_chain());
			new_hourglass.set_apex(_left.get_apaxes()[1], 0);
			r_val = concatenate_two_funnels_oc(_left.get_first_chain(), _right.get_first_chain(), common_edge, _left.get_apaxes()[0], new_hourglass.get_second_edge());
			new_hourglass.set_apex(_left.get_apaxes()[0], 1);
		}
		else {
			new_hourglass.set_first_edge(left_edge_list[0]);
			new_hourglass.set_first_chain(_left.get_first_chain());
			new_hourglass.set_apex(_left.get_apaxes()[0], 0);
			r_val = concatenate_two_funnels_oc(_left.get_second_chain(), _right.get_first_chain(), common_edge, _left.get_apaxes()[1], new_hourglass.get_second_edge());
			new_hourglass.set_apex(_left.get_apaxes()[1], 1);
		}
		if (r_val.new_string != NULL) {
			new_hourglass.set_string(new String(_left.get_string(), r_val.new_string));
			new_hourglass.set_apex(r_val.apax, 1);
		}
		new_hourglass.set_second_chain(r_val.chain);

	}
	else if (left_openess && !right_openess) {
		new_hourglass.set_string(_right.get_string());
		if (common_edge_check[0] == 0) {
			new_hourglass.set_first_edge(left_edge_list[1]);
		}
		else {
			new_hourglass.set_first_edge(left_edge_list[0]);
		}
		Return_val r_val;

		if (common_edge_check[1] == 0) {
			new_hourglass.set_second_edge(right_edge_list[1]);
			new_hourglass.set_second_chain(_right.get_second_chain());
			new_hourglass.set_apex(_right.get_apaxes()[1], 1);
			r_val = concatenate_two_funnels_oc(_right.get_first_chain(), _left.get_first_chain(), common_edge, _right.get_apaxes()[0], new_hourglass.get_first_edge());
			new_hourglass.set_apex(_right.get_apaxes()[0], 0);
		}
		else {
			new_hourglass.set_second_edge(right_edge_list[0]);
			new_hourglass.set_second_chain(_right.get_first_chain());
			new_hourglass.set_apex(_right.get_apaxes()[0], 1);
			r_val = concatenate_two_funnels_oc(_right.get_second_chain(), _left.get_first_chain(), common_edge, _right.get_apaxes()[1], new_hourglass.get_first_edge());
			new_hourglass.set_apex(_right.get_apaxes()[1], 0);
		}
		if (r_val.new_string != NULL) {
			new_hourglass.set_string(new String(r_val.new_string, _right.get_string()));
			new_hourglass.set_apex(r_val.apax, 0);
		}
		new_hourglass.set_first_chain(r_val.chain);
	}
	else {//two open hourglasses-> could be open or closed

		Hourglass result = concatenateOpenOpen(_left, _right);/////////////////////////////////////debugging�뵵�Ӵϴ�.
		return result;
		/*
		if (common_edge_check[0] == 0) {//0��° edge�� common�϶�
		new_hourglass.set_first_edge(left_edge_list[1]);//common�ƴ� �ָ� new H�� first edge�� �����
		}
		else {
		new_hourglass.set_first_edge(left_edge_list[0]);//left�� common�ƴ� edge�� newH�� first edge�� �����.
		}
		if (common_edge_check[1] == 0) {
		new_hourglass.set_second_edge(right_edge_list[1]);//right H�� common �ƴ� edge�� new H�� second edge�� �����
		}
		else {
		new_hourglass.set_second_edge(right_edge_list[0]);
		}

		//set outer chain
		Chain** left = _left.get_first_chain(), **right = _right.get_first_chain();
		int outer_chain[2][2];
		int p1 = common_edge->get_origin();

		//left[0] -> first chain of _left
		//left[1] -> second chain of _left
		//right[0] -> first chain of _right
		//right[1] -> second chain of _right

		if ((left[0]->check_inclusive(p1) && left[1]->check_inclusive(p1)) || (right[0]->check_inclusive(p1) && right[1]->check_inclusive(p1))){
		p1 = common_edge->get_dest();//left chain�� first/second chain�� point list Ȯ���ؼ� p1(common_edge�� origin)�� �Ѵ� ���Ե��ִ���, right chain�� ���ؼ��� �Ȱ��� Ȯ��
		//�� �� �ϳ��� �����ϸ�? p1 �� ���� common_edge�� destination���� ������ ��.
		//������ �ϸ�...?funnel �̶�� �� �ƴѰ�
		}
		vector<int> left_0_point_list = left[0]->get_point_list();//list of points in first chain of _left
		if (left_0_point_list.front() == p1 || left_0_point_list.back() == p1) {//if first chain of _left either start or end with p1
		outer_chain[0][0] = 0;//= left[0] -> set outer first outer chain
		outer_chain[1][0] = 1;// left[1];
		}
		else {
		outer_chain[1][0] = 0;// left[0];
		outer_chain[0][0] = 1;// left[1];
		}
		vector<int> right_0_point_list = right[0]->get_point_list();//list of points in first chain of _right
		if (right_0_point_list.front() == p1 || right_0_point_list.back() == p1) {
		outer_chain[0][1] = 0;// right[0];
		outer_chain[1][1] = 1;// right[1];
		}
		else {
		outer_chain[1][1] = 0;// right[0];
		outer_chain[0][1] = 1;// right[1];
		}

		//compute outer_chain
		Hourglass outer_hour;
		bool check = true;
		Chain *left_c_list[2];
		Chain *right_c_list[2];
		int apex[2][2];
		for (int i = 0; i < 2; i++) {
		Chain * left_upper_chain = left[outer_chain[i][0]];
		Chain * right_upper_chain = right[outer_chain[i][1]];
		Chain * left_down_chain = left[outer_chain[(i+1)%2][0]];
		Chain * right_down_chain = right[outer_chain[(i + 1) % 2][1]];
		int t1, t2;

		//���� ���⼭���� outer tangent ����ؼ� valid���� Ȯ���ϰ�, apex�� chain�� ��� index�� �ִ��� �̾Ƴ��ߵ�!////////////////////////////////////
		bool valid = compute_outer_tangent(left_upper_chain, right_upper_chain, &t1, &t2, common_edge, left_down_chain, right_down_chain);

		if (valid) {//set apex and left & right c_lists for the upcoming new chain!
		apex[i][0] = left_upper_chain->get_point(t1); //left_upper_chain�� c_point_list ���� index�� 't1'�� element�� return ��
		left_c_list[i] = left_upper_chain->cutting_chain(common_edge, apex[i][0],left_upper_chain);
		apex[i][1] = right_upper_chain->get_point(t2);
		right_c_list[i] = right_upper_chain->cutting_chain(common_edge, apex[i][1], right_upper_chain);
		}
		else {
		check = false;
		}
		}
		if (check) {
		Chain * new_chain[2];
		for (int i = 0; i < 2; i++) {
		Chain * c1 = left_c_list[i];
		Chain * c2 = right_c_list[i];
		int p1 = apex[i][0];
		int p2 = apex[i][1];
		new_chain[i] = new Chain(c1, c2, p1, p2);
		}
		new_hourglass.set_first_chain(new_chain);
		new_hourglass.duplicate_strings();
		outer_hour = new_hourglass;
		//return outer_hour;
		}
		//////////////////////////////////////////////////////


		//compute inner chain
		Hourglass h[2];
		bool h_valid[2] = { true, true };
		for (int i = 0; i < 2; i++) {

		Chain * left_upper_chain = left[outer_chain[i][0]];
		Chain * right_upper_chain = right[outer_chain[i][1]];
		Chain * left_down_chain = left[outer_chain[(i + 1) % 2][0]];
		Chain * right_down_chain = right[outer_chain[(i + 1) % 2][1]];
		int t1, t2;

		bool check = compute_inner_tangent(left_upper_chain, right_down_chain, &t1, &t2, common_edge, left_down_chain, right_upper_chain);
		if (check == false) {
		h_valid[i] = false;
		continue;
		}
		int apax1 = left_upper_chain->get_point(t1), apax2 = right_down_chain->get_point(t2);

		int apax[2] = { apax1, apax2 };
		Hourglass hh[2];
		bool hh_valid[2] = { true, true };
		for (int k = 0; k < 2; k++) {
		new_hourglass.clear_chains();
		new_hourglass.clear_apaxes();
		new_hourglass.clear_string();
		apax1 = apax[0];
		apax2 = apax[1];
		//&&side_check(left_upper_chain, right_down_chain, left_down_chain, right_upper_chain, apx)
		if (k == 0 && right_upper_chain->check_inclusive(apax1)) {
		int edge_point = right_upper_chain->get_point(0);
		if (common_edge->check_same_point(edge_point) != -1) {
		edge_point = right_upper_chain->get_last_point();
		}
		if (right_upper_chain->check_sequence(apax1, edge_point)) {
		apax2 = apax1;
		}
		}
		if (k == 1 && left_down_chain->check_inclusive(apax2)) {
		int edge_point = left_down_chain->get_point(0);
		if (common_edge->check_same_point(edge_point) != -1) {
		edge_point = left_down_chain->get_last_point();
		}
		if (left_down_chain->check_sequence(apax2, edge_point)) {
		apax1 = apax2;
		}
		}
		Chain * cutting_chain_u;
		Chain * cutting_chain_d;
		Chain * ret[2];
		cutting_chain_u = left_upper_chain->cutting_chain(common_edge, apax1, left_upper_chain);
		cutting_chain_d = left_down_chain->cutting_chain(common_edge, apax1, left_upper_chain);
		ret[0] = cutting_chain_u;
		ret[1] = cutting_chain_d;
		new_hourglass.set_first_chain(ret);
		new_hourglass.set_apex(apax1, 0);
		new_hourglass.set_string(new String(apax1, apax2));

		cutting_chain_u = right_upper_chain->cutting_chain(common_edge, apax2, right_down_chain);
		cutting_chain_d = right_down_chain->cutting_chain(common_edge, apax2, right_down_chain);
		ret[0] = cutting_chain_u;
		ret[1] = cutting_chain_d;
		new_hourglass.set_second_chain(ret);
		new_hourglass.set_apex(apax2, 1);

		if (new_hourglass.check_valid()) {
		new_hourglass.duplicate_strings();
		hh[k] = new_hourglass;
		}
		else {
		hh_valid[k] = false;
		}
		}
		point_type min_len = -1;
		for (int k = 0; k < 2; k++) {
		if (hh_valid[k] && (min_len == -1 || min_len > hh[k].get_len())) {
		h[i] = hh[k];
		min_len = hh[k].get_len();
		}
		}
		if (min_len == -1) {
		h_valid[i] = false;
		}

		/*if (right_upper_chain->check_inclusive(apax1)) {
		apax2 = apax1;
		}

		if (left_down_chain->check_inclusive(apax2)) {
		apax1 = apax2;
		}

		Chain * cutting_chain_u;
		Chain * cutting_chain_d;
		Chain * ret[2];
		cutting_chain_u = left_upper_chain->cutting_chain(common_edge, apax1, left_upper_chain);
		cutting_chain_d = left_down_chain->cutting_chain(common_edge, apax1, left_upper_chain);
		ret[0] = cutting_chain_u;
		ret[1] = cutting_chain_d;
		new_hourglass.set_first_chain(ret);
		new_hourglass.set_apax(apax1, 0);
		new_hourglass.set_string(new String(apax1, apax2));

		cutting_chain_u = right_upper_chain->cutting_chain(common_edge, apax2, right_down_chain);
		cutting_chain_d = right_down_chain->cutting_chain(common_edge, apax2, right_down_chain);
		ret[0] = cutting_chain_u;
		ret[1] = cutting_chain_d;
		new_hourglass.set_second_chain(ret);
		new_hourglass.set_apax(apax2, 1);

		if (new_hourglass.check_valid()) {
		new_hourglass.duplicate_strings();
		h[i] = new_hourglass;
		}
		else {
		h_valid[i] = false;
		}

		}
		Hourglass min;
		point_type min_len =-1;
		for (int i = 0; i < 2; i++) {
		if (h_valid[i] && (min_len == -1 || min_len > h[i].get_len())) {
		min = h[i];
		min_len = h[i].get_len();
		}
		}

		if (check && (min_len == -1 || min_len > outer_hour.get_len())) {
		min = outer_hour;
		min_len = outer_hour.get_len();
		}
		if (min_len == -1) {
		cout << "Somthings wrong" << endl;
		}
		return min;*/
	}
	return new_hourglass;

}
