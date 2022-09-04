#pragma once
#include <vector>
#include <thread>
#include <iostream>
#include <bitset>
#include <queue>
#include <array>

// code from LAP project of Stefan Guthe

template <class T>
class CacheListNode
{
public:
	T prev, next;
};

// segmented least recently used
class CacheSLRU
{
protected:
	std::vector<CacheListNode<int>> list;
	int first[2];
	int last[2];
	std::vector<int> id;
	std::vector<char> priv;
	int priv_avail;
	long long chit, cmiss;
	std::vector<int> map;

	inline void remove_entry(int i)
	{
		int l = priv[i];
		int prev = list[i].prev;
		int next = list[i].next;
		if (prev != -1) list[prev].next = next;
		if (next != -1) list[next].prev = prev;
		if (first[l] == i) first[l] = next;
		if (last[l] == i) last[l] = prev;
	}

	inline void remove_first(int i)
	{
		int l = priv[i];
		first[l] = list[i].next;
// 			if (first[l] < 0 || first[l] >= list.size()) std::cout << i << " " << first[l] << std::endl;
		list[first[l]].prev = -1;
	}

	inline void push_last(int i)
	{
		int l = priv[i];
		list[i].prev = last[l];
		list[i].next = -1;
		if (last[l] == -1)
		{
			first[l] = i;
		}
		else
		{
			list[last[l]].next = i;
		}
		last[l] = i;
	}

public:
	CacheSLRU()
	{
		first[0] = first[1] = last[0] = last[1] = -1;
		priv_avail = 0;
		chit = cmiss = 0ll;
	}
	~CacheSLRU() {}

	inline void setSize(int entries, int dim)
	{
		map.resize(dim, -1);
		list.resize(entries);
		id.resize(entries);
		priv.resize(entries);
		restart();
	}
inline void setDim(int dim)
{
	map.resize(dim, -1);
}

	inline int get_id(int idx) const
	{
		return id[idx];
	}

	inline bool find(int &idx, int i, int *oldi = nullptr)
	{
		if (map[i] == -1)
		{
			// replace
			idx = first[0];
			if (oldi != nullptr && id[idx] != -1) *oldi = id[idx];
			if (id[idx] != -1) map[id[idx]] = -1;
			id[idx] = i;
			map[i] = idx;
			remove_first(idx);
			push_last(idx);
			cmiss++;
			return false;
		}
		else
		{
			idx = map[i];
			if (priv[idx] == -1)
			{
				priv[idx] = 0;
				remove_entry(idx);
			}
			else
			{
				remove_entry(idx);
				if (priv[idx] == 0)
				{
					priv[idx] = 1;
					if (priv_avail > 0)
					{
						priv_avail--;
					}
					else
					{
						int idx1 = first[1];
						remove_first(idx1);
						priv[idx1] = 0;
						push_last(idx1);
					}
				}
			}
			push_last(idx);
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
		int entries = (int)list.size();
		int dim = (int)map.size();
// 			for (int i = 0; i < dim; i++) map[i] = -1;
		for (int i = 0; i < entries; i++) { list[i].prev = i - 1; list[i].next = (i == entries - 1) ? -1 : i + 1; }
		for (int i = 0; i < entries; i++) id[i] = -1;
		for (int i = 0; i < entries; i++) priv[i] = 0;
		chit = cmiss = 0;
		first[0] = 0;
		last[0] = entries - 1;
		first[1] = -1;
		last[1] = -1;
		priv_avail = entries >> 1;
	}

	inline void getHitMiss(long long &hit, long long &miss) { hit = chit; miss = cmiss; chit = 0; cmiss = 0; }
	inline int getEntries() { return (int)list.size(); }
}; 
