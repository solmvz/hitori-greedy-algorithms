#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

struct state
{
	vector<vector<pair<char, char>>> matrix;
	float h = 100;
	int g = 0;
	int f = 0;
	float head = 0;
	float tail = 0;
};

void build_input(state& init) //initial state
{
	ifstream input;
	input.open("sample2.txt");

	int m, n;
	input >> m; //rows
	input >> n; //columns

	init.matrix.resize(n);
	for (int i = 0; i < n; ++i)
		init.matrix[i].resize(m);

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
		{
			char c; input >> c;
			init.matrix[i][j].first = c;
			init.matrix[i][j].second = 'W';
		}
	input.close();
	return;
}
void build_state(state& S, int size, state S_origin) //builds a copy of given state
{
	S.matrix.resize(size);
	for (int i = 0; i < size; i++)
		S.matrix[i].resize(size);
	S = S_origin;
}
void print_state(state S)
{
	for (int i = 0; i < S.matrix.size(); i++)
	{
		for (int j = 0; j < S.matrix[i].size(); j++)
			cout << "(" << S.matrix[i][j].first << "," << S.matrix[i][j].second << ")";
		cout << endl;
	}
	cout << endl;
}
void output_goal_state(state S)
{
	ofstream output;
	output.open("result.txt");

	for (int i = 0; i < S.matrix.size(); i++)
	{
		for (int j = 0; j < S.matrix[i].size(); j++)
			output << S.matrix[i][j].first << " ";
		output << "\n";
	}
	output << "\n";
	for (int i = 0; i < S.matrix.size(); i++)
	{
		for (int j = 0; j < S.matrix[i].size(); j++)
			output << S.matrix[i][j].second << " ";
		output << "\n";
	}
}

bool check_for_duplicates(state S) //returns true if there are duplicates
{
	char flag_number, flag_color;
	//row
	for (int i = 0; i < S.matrix.size(); i++)
	{
		//cout << "*" << endl;
		for (int j = 0; j < S.matrix.size(); j++)
		{
			for (int r = 0; r < S.matrix.size(); r++)
			{
				if (r != i)
					if (S.matrix[i][j].first == S.matrix[r][j].first) //it is duplicate
						if(S.matrix[i][j].second!='B'&& S.matrix[r][j].second!='B') //if none of them is black
						{
							return false;
						}
			}
			for (int c = 0; c < S.matrix.size(); c++)
			{
				if (c != j)
					if (S.matrix[i][j].first==S.matrix[i][c].first)
						if (S.matrix[i][j].second != 'B' && S.matrix[i][c].second != 'B')
						{
							return false;
						}
			}
		}
	}
	return true;
}
bool check_for_adjacent_blacks(state S) //returns true if there are 2 adjacent blacks
{
	for (int i = 0; i < S.matrix.size(); i++) //two blacks beside in a row
	{
		for (int j = 0; j < S.matrix[i].size(); j++)
		{
			if (S.matrix[i][j].second == 'B')
				if (j + 1 < S.matrix.size())
					if (S.matrix[i][j + 1].second == 'B')
						return true;
		}
	}
	for (int i = 0; i < S.matrix.size(); i++) //two blacks beside in a column
	{
		for (int j = 0; j < S.matrix[i].size(); j++)
		{
			if (S.matrix[i][j].second == 'B')
				if (i + 1 < S.matrix.size())
					if (S.matrix[i + 1][j].second == 'B')
						return true;
		}
	}
	return false;
}
bool check_for_blocks(state S) //returns true if there are 4 blacks around a white
{
	int Size = S.matrix.size() - 1;
	for (int i = 0; i < S.matrix.size(); i++)
	{
		for (int j = 0; j < S.matrix.size(); j++)
		{
			if (i == 0 && j == 0) //up left corner
				if (S.matrix[i][j + 1].second == 'B' && S.matrix[i + 1][j].second == 'B')
					return true;
			if (i == 0 && j > 0 && j + 1 < S.matrix.size()) //up row
				if (S.matrix[i][j - 1].second == 'B' && S.matrix[i][j + 1].second == 'B' && S.matrix[i + 1][j].second == 'B')
					return true;
			if (i == 0 && j == Size) //up right corner
				if (S.matrix[i][j - 1].second == 'B' && S.matrix[i + 1][j].second == 'B')
					return true;
			if (i == Size && j == 0) //down left corner
				if (S.matrix[i - 1][j].second == 'B' && S.matrix[i][j + 1].second == 'B')
					return true;
			if (i == Size && j > 0 && j + 1 < S.matrix.size()) //down row
				if (S.matrix[i - 1][j].second == 'B' && S.matrix[i][j - 1].second && S.matrix[i][j + 1].second == 'B')
					return true;
			if (i == Size && j == Size) //down right corner
				if (S.matrix[i - 1][j].second == 'B' && S.matrix[i][j - 1].second == 'B')
					return true;

			if (j == 0 && i > 0 && i + 1 < S.matrix.size()) //leftmost col
				if (S.matrix[i - 1][j].second == 'B' && S.matrix[i][j + 1].second == 'B' && S.matrix[i + 1][j].second == 'B')
					return true;
			if (j == Size && i > 0 && i + 1 < S.matrix.size()) //rightmost col
				if (S.matrix[i][j - 1].second == 'B' && S.matrix[i - 1][j].second == 'B' && S.matrix[i + 1][j].second == 'B')
					return true;
			if (j > 0 && i > 0 && j + 1 < S.matrix.size() && i + 1 < S.matrix.size())
			{
				if (S.matrix[i][j - 1].second == 'B' && S.matrix[i][j + 1].second == 'B')
					if (S.matrix[i - 1][j].second == 'B' && S.matrix[i + 1][j].second == 'B')
						return true;
			}
		}
	}
	return false;
}

void color_black(state& S, int x, int y)
{
	S.matrix[x][y].second = 'B';
}
state color_state(state S, int x_black, int y_black)
{
	state new_state;
	build_state(new_state, S.matrix.size(), S);
	color_black(S, x_black, y_black);
	return S;
}

//Successor
bool is_valid(state S)
{
	if (!check_for_adjacent_blacks(S)) //state validation
		if (!check_for_blocks(S)) //false
			return true;
	return false;
}
bool search_in_successor(vector<state>succ, state S)
{
	for (int i = 0; i < succ.size(); i++)
	{
		if (succ[i].matrix == S.matrix)
			return true;
	}
	return false;
}
vector<state> successor(state S)
{
	vector<state> result;
	state new_state;

	for (int i = 0; i < S.matrix.size(); i++)
		for (int j = 0; j < S.matrix.size(); j++)
		{
			for (int r = 0; r < S.matrix.size(); r++)
            {
				if (r != i)
					if (S.matrix[i][j].first == S.matrix[r][j].first)
						if (S.matrix[i][j].second != 'B' && S.matrix[r][j].second != 'B')
						{
							new_state = color_state(S, i, j);
							if (is_valid(new_state))
							{
								if (search_in_successor(result, new_state) == false)
									result.push_back(new_state);
							}
							new_state = color_state(S, r, j);
							if (is_valid(new_state))
							{
								if (search_in_successor(result, new_state) == false)
									result.push_back(new_state);
							}
						}
			}
			for (int c = 0; c < S.matrix.size(); c++)
            {
				if (c != j)
					if (S.matrix[i][j].first == S.matrix[i][c].first)
						if (S.matrix[i][j].second != 'B' && S.matrix[i][c].second != 'B')
						{
							new_state = color_state(S, i, j);
							if (is_valid(new_state))
							{
								if (search_in_successor(result, new_state) == false)
									result.push_back(new_state);
							}
							new_state = color_state(S, i, c);
							if (is_valid(new_state))
							{
								if (search_in_successor(result, new_state) == false)
									result.push_back(new_state);
							}
						}
			}
		}
	return result;
}

//Heuristic
void init_false(vector<bool> f)
{
	for (int i = 0; i < f.size(); i++)
		f[i] = false;
}
void heuristic_count(state& S)
{
	int Count = 0;
	int c = 0;
	char flag;

	vector<bool>flags;
	flags.resize(9);
	init_false(flags);

	//row
	for (int i = 0; i < S.matrix.size(); i++)
	{
		for (int j = 0; j < S.matrix.size(); j++)
		{
			flag = S.matrix[i][j].first;
			if (flags[(flag - '0') - 1] == false)
			{
				for (int k = 0; k < S.matrix.size(); k++)
				{
					if (flag == S.matrix[i][k].first && S.matrix[i][k].second != 'B')
						c++;
				}
			}
			if (c != 1)
				Count += c;
			c = 0;
			flags[(flag - '0') - 1] = true;
		}
		for (int i = 0; i < flags.size(); i++)
			flags[i] = false;
	}

	for (int i = 0; i < flags.size(); i++)
		flags[i] = false;

	//col
	for (int i = 0; i < S.matrix.size(); i++)
	{
		for (int j = 0; j < S.matrix.size(); j++)
		{
			flag = S.matrix[j][i].first;
			if (flags[(flag - '0') - 1] == false)
			{
				for (int k = 0; k < S.matrix.size(); k++)
				{
					if (flag == S.matrix[k][i].first && S.matrix[k][i].second != 'B')
						c++;
				}
			}
			if (c != 1)
				Count += c;
			c = 0;
			flags[(flag - '0') - 1] = true;
		}
		for (int i = 0; i < flags.size(); i++)
			flags[i] = false;
	}
	S.h = Count;
}
void heuristic_whites(state& S)
{
	int count = 0;
	for (int i = 0; i < S.matrix.size(); i++)
		for (int j = 0; j < S.matrix.size(); j++) {
			if (S.matrix[i][j].second == 'W') {
				count++;
			}
		}
	S.h = count;
}

void heuristic(state& S)
{
	heuristic_count(S);
}

//Goal
bool is_goal(state S) //returns true if there are no duplicates
{
	if (!check_for_duplicates(S)) //duplicate darim
			return false;
	return true;
}

//Greedy
struct compare_h_for_greedy
{
	bool operator()(state const& s1, state const& s2)
	{
		return s1.h > s2.h;
	}
};
bool visited(state S, vector<state>& list)
{
	for (int i = 0; i < list.size(); i++)
	{
		if (list[i].matrix == S.matrix)
			return true;
	}
	list.push_back(S);
	return false;
}
state greedy(state& current_state, vector<state>& visited_list, int& max_space, int& no_of_succ)
{
	priority_queue<state, vector<state>, compare_h_for_greedy> state_q;
	vector<state>successors = successor(current_state);
	no_of_succ += successors.size();
	int q_size = 0;

	heuristic(current_state);

	for (int i = 0; i < successors.size(); i++)
	{
		heuristic(successors[i]);
		state_q.push(successors[i]);
	}
	q_size = state_q.size();
	max_space = max(q_size, max_space);

	while (!state_q.empty())
	{
		cout << "h : " << current_state.h << endl;

		if (is_goal(current_state) == true)
		{
			max_space += visited_list.size();
			return current_state;
		}

		current_state = state_q.top();
		state_q.pop();
		if (!visited(current_state, visited_list))
		{
			successors.clear();
			successors = successor(current_state);
			no_of_succ += successors.size();
			for (int i = 0; i < successors.size(); i++)
			{
				heuristic(successors[i]);
				state_q.push(successors[i]);
			}
		}
		q_size = state_q.size();
		max_space = max(q_size, max_space);
	}
}

//A Star
struct compare_f_for_astar
{
	bool operator()(state const& s1, state const& s2)
	{
		return s1.f > s2.f;
	}
};
void calculate_f(state& S)
{
	S.f = S.g + S.h;
	return;
}
state A_star(state& current_state, vector<state>& visited_list, int& max_space, int& no_of_succ)
{
	int cost = 1;
	priority_queue<state, vector<state>, compare_f_for_astar> state_q;
	vector<state>successors = successor(current_state);
	no_of_succ += successors.size();
	heuristic(current_state);

	current_state.g = 0;
	calculate_f(current_state);

	int q_size = 0;

	for (int i = 0; i < successors.size(); i++)
	{
		heuristic(successors[i]);
		successors[i].g = cost;
		calculate_f(successors[i]);
		state_q.push(successors[i]);
	}

	q_size = state_q.size();
	max_space = max(q_size, max_space);
	while (!state_q.empty())
	{
		if (is_goal(current_state) == true)
		{
			max_space += visited_list.size();
			return current_state;
		}

		current_state = state_q.top();
		state_q.pop();

		if (!visited(current_state, visited_list))
		{
			cost += 1;
			successors.clear();
			successors = successor(current_state);
			no_of_succ += successors.size();

			for (int i = 0; i < successors.size(); i++)
			{
				heuristic(successors[i]);
				successors[i].g = cost;
				calculate_f(successors[i]);
				state_q.push(successors[i]);
			}
		}
		q_size = state_q.size();
		max_space = max(q_size, max_space);
	}
}

//Hill climbing
struct compare_h_for_HL
{
	bool operator()(state const& s1, state const& s2)
	{
		return s1.h > s2.h;
	}
};
bool visited_HC(state S, vector<state>& L_seen)
{
	for (int i = 0; i < L_seen.size(); i++)
	{
		if (L_seen[i].matrix == S.matrix)
			return true;
	}
	return false;
}

// Simple
void clear(priority_queue<state, vector<state>, compare_h_for_HL> L)
{
	priority_queue<state, vector<state>, compare_h_for_HL> empty;
	swap(L, empty);
}
state simple_HC(state& init, vector<state>& L_seen, int& max_space, int& no_of_succ)
{
	priority_queue<state, vector<state>, compare_h_for_HL> L;
	heuristic(init);
	L.push(init);
	vector<state> next_states;
	state current;
	int L_size = 0;
	while (!L.empty())
	{
		current  = L.top();
		L.pop();
		if (is_goal(current))
		{
			cout << "L_seen size : " << L_seen.size() << endl;
			max_space += L_seen.size();
			return current;
		}
		clear(L);
		next_states = successor(current);
		no_of_succ += next_states.size();
		for (int i = 0; i < next_states.size(); i++) {
			heuristic(next_states[i]);
			if (!visited_HC(next_states[i], L_seen)) {
				L.push(next_states[i]);
			}
		}
		L_size = L.size();
		max_space = max(L_size, max_space);
		L_seen.push_back(current);
	}
	cout << "L_seen : " << L_seen.size() << endl;
	cout << "Sorry! Goal not found" << endl << endl;
	return current;
}

// Random
void clear_vec(vector<state>& L)
{
	vector<state> empty;
	swap(L, empty);
}
void normilize(vector<state>& L)
{
	float temp_h;
	float sum = 0;

	for (int i = 0; i < L.size(); i++)
    { //sum of state h(n)s
		sum += L[i].h;
	}
	L[0].head = 0;
	for (int i = 0; i < L.size(); i++)
	{
		temp_h = L[i].h / sum;
		L[i].h = temp_h;
		L[i].tail = temp_h + L[i].head;
		if (i + 1 < L.size()) {
			L[i + 1].head = temp_h + L[i].head;
		}
	}
}
state rand_choose(vector<state>& L)
{
	state temp_S;
	normilize(L);
	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	for (int i = 0; i < L.size(); i++)
    {
		if (r > L[i].head&& r < L[i].tail)
		{
			temp_S = L[i];
			return temp_S;
		}
	}
}
state random_HC(state& init, vector<state>& L_seen, int& max_space, int& no_of_succ)
{
	vector<state> L;
	heuristic(init);
	L.push_back(init);
	vector<state> next_states;
	state current;
	int L_size = 0;

	while (L.size()!= 0)
	{
		if (L.size() > 1)
        {
			current = rand_choose(L);
		}
		else
		{
			current = L.front();
			L.pop_back();
		}
		if (is_goal(current))
		{
			cout << "L_seen : " << L_seen.size() << endl;
			max_space += L_seen.size();
			return current;
		}
		clear_vec(L);
		next_states = successor(current);
		no_of_succ += next_states.size();
		for (int i = 0; i < next_states.size(); i++)
		{
			heuristic(next_states[i]);
			if (!visited_HC(next_states[i], L_seen)) {
				L.push_back(next_states[i]);
			}
		}
		L_size = L.size();
		max_space = max(L_size, max_space);
		L_seen.push_back(current);
	}
	cout << "L_seen : " << L_seen.size() << endl;
	cout << "Sorry! Goal not found" << endl << endl;
	return current;
}

//Simulated Annealing
state SA_rand(vector<state>& list)
{
	int r = rand() % list.size();
	return list[r];
}
void schedule(float& T,int& k,state current)
{
	float T_0 = 200;
	float alpha = 0.3;
	T = exp(-alpha * k) * T_0;
	k++;
	return;
}
bool visited_SA(state S, vector<state>& L_seen)
{
	for (int i = 0; i < L_seen.size(); i++)
	{
		if (L_seen[i].matrix == S.matrix)
			return true;
	}
	return false;
}
state simulated_annealing(state& current, vector<state>& L_seen, int& max_space, int&no_of_succ)
{
	float T = 200;
	int k = 1;
	int vec_size = 0;
	schedule(T, k, current);

	vector<state> next_states;
	state next;

	while (T > 1)
    {
		next_states = successor(current);
		no_of_succ += next_states.size();
		vec_size = next_states.size();
		max_space = max(max_space, vec_size);
		next = SA_rand(next_states);

		heuristic(next);
		heuristic(current);

		float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		int deltaE = current.h - next.h;
		float P = exp(deltaE / T);

		if (!visited_SA(next, L_seen))
        {
			if (deltaE > 0)
			{
				current = next;
			}
			else if (r < P)
			{
				current = next;
			}
		}
		L_seen.push_back(current);
		schedule(T, k, current);
	}
	cout << "steps : " << k << endl;
	if (is_goal(current))
    {
		cout << endl << "Goal" << endl;
	}
	max_space += L_seen.size();
	return current;
}

int main(void)
{
	state init_state;
	build_input(init_state);
	cout << "Init State : " << endl;
	print_state(init_state);

	srand(static_cast <unsigned> (time(0)));

	vector<state>visited_list;
	state goal_state;

	int max_space = 0;
	int no_of_succ = 0;

	auto search_start = high_resolution_clock::now();

	goal_state = greedy(init_state, visited_list, max_space, no_of_succ);

	//goal_state = A_star(init_state, visited_list, max_space, no_of_succ);

	//goal_state = simple_HC(init_state, visited_list, max_space, no_of_succ);

	//goal_state = random_HC(init_state, visited_list, max_space, no_of_succ);

	//goal_state = simulated_annealing(init_state, visited_list, max_space, no_of_succ);

	auto search_end = high_resolution_clock::now();
	auto duration_ids = duration_cast<milliseconds>(search_end - search_start);

	cout << "\nGoal State: " << endl;
	print_state(goal_state);
	output_goal_state(goal_state);
    cout << endl;

    cout << "\nSearch Time : " << duration_ids.count() << " milliseconds" << endl;
	cout << "Number of Successors: " << no_of_succ << endl;
	cout << "Maximum Saved States: " << max_space << endl;

	return 0;
}
