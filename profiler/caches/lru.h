/*
 * Copyright (C) 2022, Max von Buelow
 * Technical University of Darmstadt - Interactive Graphics Systems Group
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include <vector>
#include <thread>
#include <iostream>
#include <bitset>
#include <queue>
#include <array>

struct CacheLRU {
	struct Node {
		int prev, next;
	};
	std::vector<Node> list;
	std::vector<int> id;
	std::vector<int> map;
	int front, back, free;
 
	inline void setSize(int entries, int dim)
	{
		map.resize(dim, -1);
		list.resize(entries);
		id.resize(entries);
		restart();
	}
	inline void setDim(int dim)
	{
		map.resize(dim, -1);
	}
	inline void restart()
	{
		int entries = (int)list.size();
// 		int dim = (int)map.size();
// 		for (int i = 0; i < dim; i++) map[i] = -1;
		for (int i = 0; i < entries; i++) { id[i] = -1; list[i].prev = i - 1; list[i].next = (i == entries - 1) ? -1 : i + 1; }
		front = -1;
		back = 0;
		free = 0;
	}
	
	inline int get_id(int idx) const
	{
// 		if (idx < 0 || idx > siz
		return id[idx];
	}
	inline bool present(int x) const
	{
// 		if (x>=map.size()) std::cout << x << " " << map.size() << std::endl;
		return map[x] != -1;
	}
	// FRONT --------- BACK
	void remove_entry(int x)
	{
// 		std::cout << "re" << std::endl;
		Node &e = list[x];
// 		if (e.prev == x || e.next == x) throw std::runtime_error("invalid prevnext ptrs");
		if (e.prev != -1) list[e.prev].next = e.next;
		if (e.next != -1) list[e.next].prev = e.prev;

		if (back == x) {/*if (e.prev == -1) std::cerr << "here! " << front << " " << back << std::endl;*/ back = e.prev; }
		if (front == x) front = e.next;
// 		if(front==back) std::cout << "rem: front=back" << std::endl;

		e.prev = -1;
		e.next = free;
		free = x;
	}
	int fullc = 0;
	int push_front()
	{
// 		std::cout << "pf" << std::endl;
		int nextfree = list[free].next;

		if (front != -1) list[front].prev = free;
		list[free].next = front;
		list[free].prev = -1;

		front = free;
		free = nextfree;
// 		if(front==back) std::cout << "push: front=back " << front << " " << back << " nextfree: " << free << "@" << std::this_thread::get_id() << " " << fullc << std::endl;

		return front;
	}
   inline bool find(int &idx, int i, int *oldi = nullptr)
	{
// 		std::cout << "f" << std::endl;
		bool p = present(i);
		if (id.size() == 1) { // direct mapped
			idx = 0;
			if (!p) {
				if (oldi != nullptr/* && id[free] != -1*/) *oldi = get_id(idx);
				if (id[idx] != -1) map[id[idx]] = -1;
				id[idx] = i;
				map[i] = idx;
			}
			return p;
		}
		if (!p) {
			if (free == -1) { // full
				++fullc;
				// evict LRU
// 		if ((std::size_t)back >= 128) throw std::runtime_error("XXX " + std::to_string(back));
// 				std::cout << back << " " << id[back] << std::endl;
				map[id[back]] = -1;
				remove_entry(back);
			}
			if (oldi != nullptr/* && id[free] != -1*/) *oldi = get_id(free);
		} else {
			idx = map[i];
			if (front == idx && back == idx) return p;
			remove_entry(idx);
		}

		idx = push_front();
		id[idx] = i;
		map[i] = idx;
		return p;
	}
}; 
