#pragma once

#include <vector>
#include <random>
#include <queue>
#include <unordered_set>

#include "ray.h"
#include "aabb.h"
#include "tri_intersect.h"
#include "conn.h"

struct Vertex {
	float pos[6];
// 	float nrm[3];

	Vertex() {}
	Vertex(float x, float y, float z) : pos{ x, y, z }
	{}
	Vertex(const Vertex &vtx) : pos{ vtx.pos[0], vtx.pos[1], vtx.pos[2], vtx.pos[3], vtx.pos[4], vtx.pos[5] }
	{}
	Vertex(const volatile Vertex &vtx) : pos{ vtx.pos[0], vtx.pos[1], vtx.pos[2], vtx.pos[3], vtx.pos[4], vtx.pos[5] }
	{}

	void fill(const VtxBuf &buf, Mesh::VtxDataType type)
	{
		if (buf.has(types::PX)) {
			pos[0] = buf.get<float>(types::PX);
			pos[1] = buf.get<float>(types::PY);
			pos[2] = buf.get<float>(types::PZ);
		}
		if (buf.has(types::NX)) {
			pos[3] = buf.get<float>(types::NX);
			pos[4] = buf.get<float>(types::NY);
			pos[5] = buf.get<float>(types::NZ);
		}
	}
};

struct MyMesh : Mesh {

	bool has_normals = false;
	std::vector<Vertex> vertices;

	std::vector<AABB> aabbs;

	void shuffle()
	{
		std::cout << "DO SHUFFLE" << std::endl;
		//bfs(); // bring faces in order
		std::random_device rd;
		std::mt19937 g(rd());

// 		std::vector<uint32_t> fperm(faces.size());
// 		std::iota(fperm.begin(), fperm.end(), 0);
// 		std::shuffle(fperm.begin(), fperm.end(), g);
// 		std::vector<Face> ftmp = faces;
// 		for (uint32_t i = 0; i < faces.size(); ++i) {
// 			faces[fperm[i]] = ftmp[i];
// 		}

		std::vector<uint32_t> vperm(vertices.size());
		std::iota(vperm.begin(), vperm.end(), 0);
		std::shuffle(vperm.begin(), vperm.end(), g);

		std::vector<Vertex> vtmp = vertices;
		for (uint32_t i = 0; i < vertices.size(); ++i) {
			vertices[vperm[i]] = vtmp[i];
		}
		for (uint32_t i = 0; i < faces.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				faces[i].idx[j] = vperm[faces[i].idx[j]];
			}
		}
	}
	void bfs()
	{
		std::cout << "DO BFS" << std::endl;
		// traverse mesh
		std::vector<orgedge> edges(faces.size() * 3);
		Conn conn(edges.data(), edges.size());
		ConnBuilder builder((const uint32_t*)faces.data(), faces.size(), conn);
		builder.build();

// 		std::cout << faces.size() << std::endl;
// 		faces.clear(); // already stored in edges
		std::vector<uint32_t> vperm(vertices.size(), -1u);
		uint32_t vidx = 0;

		std::unordered_set<uint32_t> outstanding;
		for (uint32_t i = 0; i < edges.size() / 3 /* faces.size() is 0 now */; ++i) {
			outstanding.insert(i);
		}
		while (!outstanding.empty()) {
			uint32_t f = *outstanding.begin();
// 			std::cout << "proc " << f << std::endl;
			std::queue<uint32_t> q;
			std::unordered_set<uint32_t> seen;
			q.push(f);
			do {
				uint32_t ff = q.front(); q.pop();
				if (seen.find(ff) != seen.end()) continue;
				uint32_t idx[3];
				seen.insert(ff);
				outstanding.erase(ff);
				for (int i = 0; i < 3; ++i) {
					uint32_t h = conn.mkhe(ff, i);
					uint32_t vi = conn.org(h);
					if (vperm[vi] == -1u) {
						vi = (vperm[vi] = ++vidx);
					} else {
						vi = vperm[vi];
					}
					idx[i] = vi;
					uint32_t t = conn.twin(h);
					if (t == f) continue;
					uint32_t fff = conn.face(t);
					if (seen.find(fff) != seen.end()) continue;
					q.push(fff);
				}
// 				faces.push_back(Face(idx[0], idx[1], idx[2]));
			} while (!q.empty());
		}
// 		std::cout << faces.size() << std::endl;
		std::vector<Vertex> vtmp = vertices;
		for (uint32_t i = 0; i < vertices.size(); ++i) {
			if (vperm[i] == -1u) continue;
			vertices[vperm[i]] = vtmp[i];
		}
		for (uint32_t i = 0; i < faces.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				faces[i].idx[j] = vperm[faces[i].idx[j]];
			}
		}
	}

	void precompute_aabbs()
	{
		aabbs.resize(faces.size());
		for (int i = 0; i < aabbs.size(); ++i) {
			for (int v = 0; v < 3; ++v) {
				aabbs[i].feed(vertices[faces[i].idx[v]].pos);
			}
		}
	}
	void compute_normals()
	{
		if (has_normals) return;
		compute_normals(0, size());
	}
	void compute_normals(uint32_t o, uint32_t n)
	{
		std::cout << "Computing normals " << o << " " << n << std::endl;
		for (int i = o; i < o + n; ++i) {
			float ux = get_coord(i, 1)[0] - get_coord(i, 0)[0], uy = get_coord(i, 1)[1] - get_coord(i, 0)[1], uz = get_coord(i, 1)[2] - get_coord(i, 0)[2];
			float vx = get_coord(i, 2)[0] - get_coord(i, 0)[0], vy = get_coord(i, 2)[1] - get_coord(i, 0)[1], vz = get_coord(i, 2)[2] - get_coord(i, 0)[2];

			float nx = -uy * vz + uz * vy;
			float ny = -uz * vx + ux * vz;
			float nz = -ux * vy + uy * vx;
			float nl = std::sqrt(nx * nx + ny * ny + nz * nz);
			if (nl != 0.f) { nx /= nl; ny /= nl; nz /= nl; }

			for (int j = 0; j < 3; ++j) {
				vertices[faces[i].idx[j]].pos[3] = nx;
				vertices[faces[i].idx[j]].pos[4] = ny;
				vertices[faces[i].idx[j]].pos[5] = nz;
			}
		}
	}
	AABB compute_aabb() const
	{
		AABB aabb;
		for (int i = 0; i < vertices.size(); ++i) {
			aabb.feed(vertices[i].pos);
		}
		return aabb;
	}

	uint32_t idx_counter[3] = { 0, 0, 0 };
	void add_vtx(const VtxBuf &buf, Mesh::VtxDataType type)
	{
		if (vertices.size() <= idx_counter[type]) {
			vertices.emplace_back();
		}
		if (buf.has(types::NX)) has_normals = true;
		vertices[idx_counter[type]++].fill(buf, type);
	}

	void validate(const VtxBuf &buf) const
	{
		if (!buf.has(types::PX) || !buf.has(types::PY) || !buf.has(types::PZ)) throw std::runtime_error("Invalid vertex defintion: no 3D position");
	}

	std::size_t size_vertices() const
	{
		return vertices.size();
	}
	std::size_t size() const
	{
		return faces.size();
	}
	float get_bounding_min(int axis, uint32_t idx) const
	{
		return aabbs[idx].min[axis];
	}
	float get_bounding_max(int axis, uint32_t idx) const
	{
		return aabbs[idx].max[axis];
	}

	const float *get_coord(int face, int edge) const
	{
		return vertices[faces[face].idx[edge]].pos;
	}
	bool intersect(float &t, void *data, uint32_t idx, const Ray &ray) const
	{
		float u, v;
		if (!tri_intersect(t, u, v, ray, get_coord(idx, 0), get_coord(idx, 1), get_coord(idx, 2))) return false;

		float *d = (float*)data;
		for (int i = 0; i < 6; ++i) { // NOTE: n comes after v in memory
			d[i] = vertices[faces[idx].idx[0]].pos[i] * (1.f - u - v) + vertices[faces[idx].idx[1]].pos[i] * u + vertices[faces[idx].idx[2]].pos[i] * v;
		}

		return true;
	}
};
