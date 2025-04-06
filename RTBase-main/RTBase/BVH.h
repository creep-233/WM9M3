//#pragma once
//
//#include "Core.h"
//#include "Geometry.h"
//#include <vector>
//#include <algorithm>
//#include <limits>
//
//class BVHNode {
//public:
//	AABB bounds;
//	BVHNode* left = nullptr;
//	BVHNode* right = nullptr;
//	int startIndex = 0;
//	int count = 0;
//	bool isLeaf = false;
//
//	void build(std::vector<Triangle>& triangles, int start, int end, int depth = 0) {
//		startIndex = start;
//		count = end - start;
//
//		// 计算当前节点包围盒
//		bounds = AABB();
//		for (int i = start; i < end; ++i) {
//			for (int j = 0; j < 3; ++j)
//				bounds.extend(triangles[i].vertices[j].p);
//		}
//
//		// 终止条件
//		if (count <= 4 || depth > 16) {
//			isLeaf = true;
//			return;
//		}
//
//		// 找出最长轴
//		Vec3 extent = bounds.max - bounds.min;
//		int axis = 0;
//		if (extent.y > extent.x) axis = 1;
//		if (extent.z > extent[axis]) axis = 2;
//
//		// 根据中心点在最长轴方向上排序
//		std::sort(triangles.begin() + start, triangles.begin() + end, [axis](const Triangle& a, const Triangle& b) {
//			float ca = (a.vertices[0].p[axis] + a.vertices[1].p[axis] + a.vertices[2].p[axis]) / 3.0f;
//			float cb = (b.vertices[0].p[axis] + b.vertices[1].p[axis] + b.vertices[2].p[axis]) / 3.0f;
//			return ca < cb;
//			});
//
//		int mid = (start + end) / 2;
//		left = new BVHNode();
//		right = new BVHNode();
//		left->build(triangles, start, mid, depth + 1);
//		right->build(triangles, mid, end, depth + 1);
//	}
//
//	bool intersect(const Ray& ray, float& tNear, int& triIndex, std::vector<Triangle>& triangles) {
//		float t0, t1;
//		if (!bounds.rayIntersect(ray, t0, t1)) return false;
//
//		bool hit = false;
//		if (isLeaf) {
//			for (int i = startIndex; i < startIndex + count; ++i) {
//				float t, u, v;
//				if (triangles[i].rayIntersect(ray, t, u, v)) {
//					if (t < tNear) {
//						tNear = t;
//						triIndex = i;
//						hit = true;
//					}
//				}
//			}
//		}
//		else {
//			bool hitL = left && left->intersect(ray, tNear, triIndex, triangles);
//			bool hitR = right && right->intersect(ray, tNear, triIndex, triangles);
//			hit = hitL || hitR;
//		}
//		return hit;
//	}
//
//	bool traverseVisible(const Ray& ray, const std::vector<Triangle>& triangles, float maxT) {
//		float t0, t1;
//		if (!bounds.rayIntersect(ray, t0, t1) || t0 > maxT) return false;
//
//		if (isLeaf) {
//			for (int i = startIndex; i < startIndex + count; ++i) {
//				float t, u, v;
//				if (triangles[i].rayIntersect(ray, t, u, v)) {
//					if (t < maxT) return false; // 被遮挡
//				}
//			}
//			return true; // 没被遮挡
//		}
//
//		bool leftVisible = left ? left->traverseVisible(ray, triangles, maxT) : true;
//		if (!leftVisible) return false;
//		return right ? right->traverseVisible(ray, triangles, maxT) : true;
//	}
//};


#pragma once

#include "Geometry.h"
#include <vector>
#include <algorithm>
#include <cfloat>

class BVHNode
{
public:
    AABB bounds;
    BVHNode* left = nullptr;
    BVHNode* right = nullptr;
    std::vector<int> triangleIndices;

    // build BVH
    void build(std::vector<Triangle>& triangles, const std::vector<int>& triIndices, int depth = 0)
    {
        triangleIndices = triIndices;
        bounds.reset();
        for (int i : triangleIndices)
        {
            bounds.extend(triangles[i].vertices[0].p);
            bounds.extend(triangles[i].vertices[1].p);
            bounds.extend(triangles[i].vertices[2].p);
        }

        // Termination condition
        if (triangleIndices.size() <= 4 || depth > 16)
        {
            return;
        }

        // Calculate the bounding box center and select the longest axis
        AABB centroidBounds;
        for (int i : triangleIndices)
        {
            Vec3 centroid = triangles[i].centre();
            centroidBounds.extend(centroid);
        }
        int axis = 0;
        Vec3 extent = centroidBounds.max - centroidBounds.min;
        if (extent.y > extent.x && extent.y > extent.z) axis = 1;
        else if (extent.z > extent.x) axis = 2;

        // Sort by center point
        std::vector<int> sorted = triangleIndices;
        std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
            return triangles[a].centre()[axis] < triangles[b].centre()[axis];
            });


        int mid = sorted.size() / 2;
        std::vector<int> leftIndices(sorted.begin(), sorted.begin() + mid);
        std::vector<int> rightIndices(sorted.begin() + mid, sorted.end());

        left = new BVHNode();
        right = new BVHNode();
        left->build(triangles, leftIndices, depth + 1);
        right->build(triangles, rightIndices, depth + 1);


        triangleIndices.clear();
    }


    void traverse(const Ray& ray, const std::vector<Triangle>& triangles, IntersectionData& closestHit) const
    {
        float tMin, tMax;
        if (!bounds.rayIntersect(ray, tMin, tMax)) return;

        if (left == nullptr && right == nullptr)
        {
            for (int i : triangleIndices)
            {
                float t, u, v;
                if (triangles[i].rayIntersect(ray, t, u, v) && t < closestHit.t)
                {
                    closestHit.t = t;
                    closestHit.ID = i;
                    closestHit.alpha = u;
                    closestHit.beta = v;
                    closestHit.gamma = 1.0f - u - v;
                }
            }
            return;
        }

        if (left) left->traverse(ray, triangles, closestHit);
        if (right) right->traverse(ray, triangles, closestHit);
    }

    // shadow ray
    bool traverseVisible(const Ray& ray, const std::vector<Triangle>& triangles, float maxT) const
    {
        float tMin, tMax;
        if (!bounds.rayIntersect(ray, tMin, tMax) || tMin > maxT) return false;

        if (left == nullptr && right == nullptr)
        {
            for (int i : triangleIndices)
            {
                float t, u, v;
                if (triangles[i].rayIntersect(ray, t, u, v) && t < maxT)
                    return true;
            }
            return false;
        }

        bool hitLeft = left ? left->traverseVisible(ray, triangles, maxT) : false;
        if (hitLeft) return true;
        bool hitRight = right ? right->traverseVisible(ray, triangles, maxT) : false;
        return hitRight;
    }
};

