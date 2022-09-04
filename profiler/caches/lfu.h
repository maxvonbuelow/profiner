#pragma once
#include <vector>
#include <thread>
#include <iostream>
#include <bitset>
#include <queue>
#include <array>

// code from LAP project of Stefan Guthe

// least frequently used
class CacheLFU
{
protected:
	long long chit, cmiss;
	std::vector<int> map;
	std::vector<int> id;
	std::vector<int> count;
	std::vector<int> order;
	std::vector<int> pos;
	int entries;

	inline void advance(int start)
	{
		// this uses a heap now
		bool done = false;
		int i = start;
		int ii = id[order[i]];
		int ci = count[ii];
		while (!done)
		{
			int l = i + i + 1;
			int r = l + 1;
			if (l >= entries) done = true;
			else
			{
				if (r >= entries)
				{
					int il = id[order[l]];
					int cl;
					if (il == -1) cl = -1; else cl = count[il];
					if ((ci > cl) || ((ci == cl) && (ii < il)))
					{
						std::swap(order[i], order[l]);
						pos[order[i]] = i;
						pos[order[l]] = l;
						i = l;
						ii = id[order[i]];
						ci = count[ii];
					}
					else
					{
						done = true;
					}
				}
				else
				{
					int il = id[order[l]];
					int ir = id[order[r]];
					int cl, cr;
					if (il == -1) cl = -1; else cl = count[il];
					if (ir == -1) cr = -1; else cr = count[ir];
					if ((cr > cl) || ((cr == cl) && (ir < il)))
					{
						// left
						if ((ci > cl) || ((ci == cl) && (ii < il)))
						{
							std::swap(order[i], order[l]);
							pos[order[i]] = i;
							pos[order[l]] = l;
							i = l;
						}
						else
						{
							done = true;
						}
					}
					else
					{
						// right
						if ((ci > cr) || ((ci == cr) && (ii < ir)))
						{
							std::swap(order[i], order[r]);
							pos[order[i]] = i;
							pos[order[r]] = r;
							i = r;
						}
						else
						{
							done = true;
						}
					}
				}
			}
		}
	}

public:
	CacheLFU()
	{
		entries = 0;
		chit = cmiss = 0ll;
	}
	~CacheLFU() {}

	inline int get_id(int idx) const
	{
		return id[idx];
	}
	inline void setSize(int p_entries, int dim)
	{
		entries = p_entries;
		map.resize(dim);
		count.resize(dim);
		id.resize(entries);
		order.resize(entries);
		pos.resize(entries);
		for (int i = 0; i < dim; i++) map[i] = -1;
		for (int i = 0; i < dim; i++) count[i] = 0;
		for (int i = 0; i < entries; i++) order[i] = i;
		for (int i = 0; i < entries; i++) pos[i] = i;
		for (int i = 0; i < entries; i++) id[i] = -1;
		chit = cmiss = 0;
	}

	inline bool find(int &idx, int i, int *oldi = nullptr)
	{
		if (map[i] == -1)
		{
			// replace
			idx = order[0];
			if (oldi != nullptr && id[idx] != -1) *oldi = id[idx];
			if (id[idx] != -1) map[id[idx]] = -1;
			id[idx] = i;
			map[i] = idx;
			count[i]++;
			advance(0);
			cmiss++;
			return false;
		}
		else
		{
			idx = map[i];
			count[i]++;
			advance(pos[idx]);
			chit++;
			return true;
		}
	}

	inline bool present(int i)
	{
		return (map[i] >= 0);
	}

	inline void restart()
	{
		int dim = (int)count.size();
		for (int i = 0; i < dim; i++) count[i] = 0;
	}

	inline void getHitMiss(long long &hit, long long &miss) { hit = chit; miss = cmiss; chit = 0; cmiss = 0; }
	inline int getEntries() { return entries; }
}; 
