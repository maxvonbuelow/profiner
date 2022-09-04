#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "types.h"
#include <iostream>

struct orgedge {
	uint32_t twin;
	uint32_t org;
};
struct AbstrConn {
	virtual uint32_t size_faces() const = 0;
	virtual uint32_t size_edges() const = 0;

	virtual uint32_t org(uint32_t he) const = 0;
	virtual uint32_t org(uint32_t face, uint32_t ledge) const = 0;
	virtual uint32_t face(uint32_t he) const = 0;
	virtual uint32_t ledge(uint32_t he) const = 0;
	virtual uint32_t next(uint32_t he) const = 0;
	virtual uint32_t prev(uint32_t he) const = 0;
	virtual uint32_t twin(uint32_t he) const = 0;
	virtual uint32_t mkhe(uint32_t face, uint32_t ledge = 0) const = 0;
	virtual uint32_t first_border_he() const = 0;
	virtual const std::unordered_set<uint32_t> &get_faces() const = 0;
};

struct Conn : AbstrConn {
	Conn()
	{}

	Conn(orgedge *e, uint32_t n) : edges(e), ne(n)
	{}

// 	std::vector<orgedge> edges;
	orgedge *edges;
	uint32_t ne;

	uint32_t size_faces() const
	{
		return size_edges() / 3;
	}
	uint32_t size_edges() const
	{
		return ne;
	}

	uint32_t org(uint32_t he) const
	{
		return edges[he].org;
	}
	uint32_t org(uint32_t face, uint32_t ledge) const
	{
		return org(face * 3 + ledge);
	}
	uint32_t face(uint32_t he) const
	{
		return he / 3;
	}
	uint32_t ledge(uint32_t he) const
	{
		return he % 3;
	}
	uint32_t next(uint32_t he) const
	{
		uint32_t f = face(he);
		uint32_t l = ledge(he);
		return mkhe(f, l == 2 ? 0 : l + 1);
	}
	uint32_t prev(uint32_t he) const
	{
		uint32_t f = face(he);
		uint32_t l = ledge(he);
		return mkhe(f, l == 0 ? 2 : l - 1);
	}
	uint32_t twin(uint32_t he) const
	{
		return edges[he].twin;
	}
	uint32_t mkhe(uint32_t face, uint32_t ledge = 0) const
	{
		return face * 3 + ledge;
	}
	uint32_t first_border_he() const
	{
		for (int i = 0; i < size_edges(); ++i) {
			if (twin(i) == i) return i;
		}
		return 0;
	}
	std::unordered_set<uint32_t> x;
	const std::unordered_set<uint32_t> &get_faces() const
	{
		return x;
	}
};
struct ConnBuilder {
	struct pairhash {
		template <typename T, typename U>
		std::size_t operator()(const std::pair<T, U> &x) const
		{
			return std::hash<T>()(x.first) + std::hash<U>()(x.second);
		}
	};
	typedef std::pair<uint32_t, uint32_t> edgemap_e;
	typedef std::unordered_map<edgemap_e, uint32_t, pairhash> edgemap;
	edgemap em;
	const uint32_t *faces;
	uint32_t nf;
	Conn &conn;

	ConnBuilder(const uint32_t *f, uint32_t nf, Conn &c) : faces(f), nf(nf), conn(c)
	{}

	inline void add_edge(uint32_t eid, uint32_t a, uint32_t b)
	{
// 		if (a == 21216) std::cout << "21216!!" << std::endl;
// 		if (a == b) return;
		uint32_t twinoff = eid;
		edgemap::iterator twin = em.find(edgemap_e(b, a));
		if (twin != em.end()) { // found a twin, merge
			std::swap(twinoff, conn.edges[twin->second].twin);
			em.erase(twin);
		} else {
			em.insert(std::make_pair(edgemap_e(a, b), twinoff));
		}
		conn.edges[eid] = orgedge{ twinoff, a };
	}
	inline void build()
	{
// 		conn.nf = mesh.faces.size();
		for (int i = 0; i < nf; ++i) {
			add_edge(i * 3, faces[i * 3 + 0], faces[i * 3 + 1]);
			add_edge(i * 3 + 1, faces[i * 3 + 1], faces[i * 3 + 2]);
			add_edge(i * 3 + 2, faces[i * 3 + 2], faces[i * 3 + 0]);
		}
	}
};
