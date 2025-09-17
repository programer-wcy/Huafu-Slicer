#include "../ClipperUtils.hpp"
#include "../ExPolygon.hpp"
#include "../Surface.hpp"
#include "../ShortestPath.hpp"
#include <windows.h>
#include <sstream>
#include <string>
#include <filesystem>
#include "FillBridge.hpp"
#include <algorithm>

namespace Slic3r {


    Polylines FillBridge::fill_surface(const Surface* surface, const FillParams& params)
    {
        Polylines polylines_out;
        polylines_out.clear();
        double area = surface->expolygon.contour.area();
        if (area < 900000000) {
            return polylines_out;
        }
        Points points = surface->expolygon.contour.points;
        Point first = points[0];
        points.push_back(first);


        std::vector<Point_t> outerBoundary;
        outerBoundary.reserve(points.size() + 1);
        for (Point one : points) {
            double x = one.x();
            double y = one.y();
            Point_t temp = { x,y };
            outerBoundary.push_back(temp);
        }
        bg::append(b_polygon.outer(), outerBoundary);

        Polygons inPolygons = surface->expolygon.holes;
        for (Polygon one : inPolygons) {
            Ring innerRing;
            innerRing.reserve(one.points.size() + 1);
            Points points = one.points;
            Point first = points[0];
            points.push_back(first);
            for (Point two : points) {
                double x = two.x();
                double y = two.y();
                Point_t temp = { x,y };
                //BOOST_LOG_TRIVIAL(info) << "Point:( " << x << " ," << y << ")";
                innerRing.push_back(temp);
            }
            b_polygon.inners().push_back(innerRing);
        }
        try {
            generateRings();
            formatTree();
            std::vector<IdIndex> visited;
            dfs(findNode({ 1,1 }), visited);
            std::reverse(bridges.begin(), bridges.end());

            //形成图
            for (auto& id : visited) {
                RingNode& node = findNode(id); //环
                std::vector<std::pair<Point_t, Point_t>> bridge_pairs; //环上桥接点对
                std::vector<Point_t> ring_points;
                std::vector<Vertex> ring_verteies;
                std::vector<int> _point_index_list;
                remove_duplicate_points(node.ring); //去掉环中重复的点

                //本环上的桥接点对
                for (auto& bridge : bridges) {
                    if (bridge.from_ii == id) {
                        bridge_pairs.push_back(std::make_pair(bridge.from, bridge.from2));
                    }
                    if (bridge.to_ii == id) {
                        bridge_pairs.push_back(std::make_pair(bridge.to, bridge.to2));
                    }
                }
 
                //需删除桥接点对的点
                for (auto& pair : bridge_pairs) {
                    int h, k;
                    for (int i=0; i < node.ring.size(); i++) {
                        if (bg::equals(node.ring[i], pair.first)) {
                            k = i;
                        }else if (bg::equals(node.ring[i], pair.second)) {
                            h = i;
                        }
                    }
                    int start = k < h ? k : h;
                    int end = k > h ? k : h;
                    if (end - start < node.ring.size() - end + start) {
                        for (int j = start + 1; j < end; j++) {
                            _point_index_list.push_back(j);
                        }
                    }
                    else {
                        for (int j = 0; j < start; j++) {
                            _point_index_list.push_back(j);
                        }
                        for (int j = end + 1; j < node.ring.size(); j++) {
                            _point_index_list.push_back(j);
                        }
                    }
                }
                for (int i = 0; i < node.ring.size(); i++) {
                    if (!isInVertexList(i, _point_index_list)) {
                        ring_points.push_back(node.ring[i]);
                    }
                }
                //形成顶点
                for (int i = 0; i < ring_points.size(); i++) {
                    Vertex v = boost::add_vertex(ring_points[i], graph);
                    ring_verteies.push_back(v);
                }
                
                //形成边
                for (int i = 0; i < ring_verteies.size() - 1; i++) {
                    boost::add_edge(ring_verteies[i], ring_verteies[i+1], graph);
                }
                boost::add_edge(ring_verteies[ring_verteies.size() - 1], ring_verteies[0], graph);

                //删除桥接点对的边
                for (auto& pair : bridge_pairs) {
                    Vertex v_s, v_e;
                    bool s = false, e = false;
                    for (Vertex v : ring_verteies) {
                        if (bg::equals(graph[v], pair.first)) {
                            v_s = v;
                            s = true;
                        }
                        else if (bg::equals(graph[v], pair.second)) {
                            v_e = v;
                            e = true;
                        }
                        if (s && e) {
                            break;
                        }
                    }
                    // 检查边是否存在
                    std::pair<Edge, bool> edge_pair = boost::edge(v_s, v_e, graph);
                    if (edge_pair.second) {
                        // 边存在，删除它
                        boost::remove_edge(edge_pair.first, graph);
                    }
                }
                RingVertices rv(id, bridge_pairs, ring_verteies);
                rv_vertices.push_back(rv);
            }

            //形成桥接边
            for (auto& bridge : bridges) {
                Vertex v1, v2, v3, v4;
                bool from = false, to = false;
                for (auto& rv : rv_vertices) {
                    if (rv.ii == bridge.from_ii) {
                        v1 = find_vertex_by_point(rv.ring_vertices, graph, bridge.from);
                        v3 = find_vertex_by_point(rv.ring_vertices, graph, bridge.from2);
                        from = true;
                    }
                    if (rv.ii == bridge.to_ii) {
                        v2 = find_vertex_by_point(rv.ring_vertices, graph, bridge.to);
                        v4 = find_vertex_by_point(rv.ring_vertices, graph, bridge.to2);
                        to = true;
                    }
                    if (from && to) {
                        break;
                    }
                }
                if (from && to) {
                    //桥接边
                    boost::add_edge(v1, v2, graph); //构成边
                    boost::add_edge(v3, v4, graph); //构成边
                }

            }
            removeDegreeOneVertices(graph);
            Vertex start = rv_vertices[0].ring_vertices[0];
            std::vector<Vertex> vertices = findEulerPath(graph, start);

            size_t s = 0;
            for (; s < vertices.size() - 1; s++) {
                Point_t p1 = graph[vertices[s]];  // 起点坐标
                Point_t p2 = graph[vertices[s + 1]];  // 终点坐标
                b_polylines.push_back({ { p1.x(),p1.y() }, { p2.x() ,p2.y() } });
            }
            Point_t p1 = graph[vertices[vertices.size() -1]];  // 起点坐标
            Point_t p2 = graph[vertices[0]];  // 终点坐标
            b_polylines.push_back({ { p1.x(),p1.y() }, { p2.x() ,p2.y() } });
        }
        catch (const std::exception& e) {
            std::string error = std::string(e.what());
        }
        return b_polylines; // polylines_out; //
    }


    //判断两个点是否相等
    bool FillBridge::equal(Point_t p1, Point_t p2) {
        return p1.x() == p2.x() && p1.y() == p2.y();
    }

    //将节点加入处理节点集合中
    RingNode FillBridge::formatNode(const size_t id, const size_t index, const Ring& ring, const int orientation) {
        Ring ringNode;
        ringNode.assign(ring.begin(), ring.end());
        RingNode node({ id, index }, ringNode, orientation);
        return node;
    }
    void FillBridge::addNode(RingNode& node) {

        // 检查键是否存在
        const auto it = ringNodes.find(node.id.id);
        if (it != ringNodes.end()) {//键就是Exy中的x，ringnodes是一个map，前面就应该是x，后面是所有x的环
            // 键已存在，将新节点添加到对应向量的末尾
            it->second.emplace_back(node);
        }
        else {
            // 键不存在，创建一个新的向量并插入到 map 中
            ringNodes[node.id.id] = { node };//ringnodes看来是最终的环集
        }
    }

    void FillBridge::removeNode(IdIndex id) {
        auto it = ringNodes.find(id.id);
        if (it != ringNodes.end()) {
            std::vector<RingNode>& nodes = it->second;
            nodes.erase(std::remove_if(nodes.begin(), nodes.end(),
                [id](RingNode x) { return x.id.id == id.id && x.id.index == id.index; }),
                nodes.end());
        }
    }

    //查找节点
    RingNode& FillBridge::findNode(const IdIndex ii) {
        const auto it = ringNodes.find(ii.id);
        std::vector<RingNode>& nodes = it->second;
        const auto it1 = std::find_if(nodes.begin(), nodes.end(), [ii](const RingNode& node) {
            return node.id.index == ii.index;
            });
        return *it1;
    }


    // 判断多边形A是否包含多边形B
    bool FillBridge::isPolygonContained(const Polygon_t& polyA, const Polygon_t& polyB) {
        // 检查B的所有点是否都在A内部
        for (const auto& point : bg::exterior_ring(polyB)) {
            if (!bg::within(point, polyA)) {
                return false;
            }
        }

        // 检查A和B的环是否有交叉
        std::vector<Polygon_t> intersection;
        bg::intersection(polyA, polyB, intersection);

        // 如果交集的总面积等于B的面积，则A包含B
        double totalArea = 0.0;
        for (const auto& poly : intersection) {
            totalArea += bg::area(poly);
        }

        return std::abs(totalArea - bg::area(polyB)) < 1e-9;
    }

    bool FillBridge::polyIntersect(const Polygon_t poly1, const Polygon_t poly2)
    {
        // 首先检查是否相交（包括边界相交）
        if (!bg::intersects(poly1, poly2)) {
            return false;
        }

        if (bg::area(poly1) == bg::area(poly2)) {
            return false;
        }
        // 如果存在包含关系，则不算作相交
        if (isPolygonContained(poly1, poly2) || isPolygonContained(poly2, poly1)) {
            return false;
        }

        return true;
    }


    //获取环偏移后的环（可能为空）
    std::vector<Ring> FillBridge::offsetRing(const Ring& ring, const double distance, double area_threshold) {
        // 将环包装成多边形
        Polygon_t inputPolygon;
        inputPolygon.outer() = ring;
        // 修复几何数据
        bg::correct(inputPolygon);

        // 使用 multi_polygon 作为 buffer 的输出类型
        bg::model::multi_polygon<Polygon_t> result;
        // 使用对称距离策略
        bg::buffer(inputPolygon, result,
            bg::strategy::buffer::distance_symmetric<double>(distance),
            bg::strategy::buffer::side_straight(),
            bg::strategy::buffer::join_round(),
            bg::strategy::buffer::end_round(),
            bg::strategy::buffer::point_circle());
        // 提取偏移后的环
        std::vector<Ring> offsetRing;
        for (const auto& poly : result) {
            if (std::abs(bg::area(poly.outer())) > area_threshold) {
                offsetRing.emplace_back(poly.outer());
            }
        }
        return offsetRing;
    }

    auto FillBridge::get_bbox(const Polygon_t& poly) {
        bg::model::box<Point_t> bbox;
        bg::envelope(poly, bbox);
        return bbox;
    }

    // 改进的安全合并函数
    Polygon_t FillBridge::safe_union(const Polygon_t& a, const Polygon_t& b) {
        const double tolerance = 1e-6;

        try {
            // 确保多边形有效
            Polygon_t valid_a = a;
            Polygon_t valid_b = b;

            if (!bg::is_valid(valid_a)) {
                bg::correct(valid_a);
            }

            if (!bg::is_valid(valid_b)) {
                bg::correct(valid_b);
            }

            // 首先检查两个多边形是否相交
            if (!polyIntersect(valid_a, valid_b)) {
                // 如果不相交，尝试合并它们
                MultiPolygon union_result;
                bg::union_(valid_a, valid_b, union_result);

                // 如果合并成功，返回结果
                if (union_result.size() == 1) {
                    return union_result[0];
                }

                // 如果合并产生多个多边形，尝试使用凸包
                if (union_result.size() > 1) {
                    Polygon_t hull;
                    bg::convex_hull(union_result, hull);

                    if (bg::is_valid(hull)) {
                        return hull;
                    }
                }

                // 如果所有方法都失败，返回第一个多边形
                return valid_a;
            }

            // 如果多边形相交，尝试合并它们
            MultiPolygon union_result;
            bg::union_(valid_a, valid_b, union_result);

            // 如果合并成功，返回结果
            if (union_result.size() == 1) {
                return union_result[0];
            }

            // 如果合并产生多个多边形，尝试使用交集和差异来合并
            if (union_result.size() > 1) {
                // 尝试找到最大的多边形
                double max_area = 0;
                Polygon_t largest_poly;

                for (const auto& poly : union_result) {
                    double area = bg::area(poly);
                    if (area > max_area) {
                        max_area = area;
                        largest_poly = poly;
                    }
                }

                return largest_poly;
            }

            // 如果所有方法都失败，返回第一个多边形
            return valid_a;
        }
        catch (const std::exception& e) {
            std::cerr << "Precise union failed: " << e.what() << std::endl;
            return a;
        }
    }

    // 改进的多边形合并函数
    void FillBridge::polygonsMerge(
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> inputs,
        std::vector<std::pair<Polygon_t, IdIndex>>& outputs, size_t orientation) {

        // 使用R-tree加速多边形合并
        typedef std::pair<bg::model::box<Point_t>, size_t> BoxValue;
        bgi::rtree<BoxValue, bgi::quadratic<16>> rtree;

        for (size_t i = 0; i < inputs.size(); i++) {
            rtree.insert(std::make_pair(get_bbox(inputs[i].first), i));
        }//构建rtree

        std::vector<bool> merged(inputs.size(), false);//初始化标识全部是未被处理
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> combined_inner_polys;

        // 首先，尝试合并所有明显相交的多边形
        for (size_t i = 0; i < inputs.size(); i++) {
            if (merged[i]) continue;//检查该环是否已经被处理

            auto current = inputs[i];
            merged[i] = true;//标记该环已被处理

            // 查找所有可能相交的多边形
            std::vector<BoxValue> candidates;//储存查找到的与本多边形相交的多边形
            auto bbox = get_bbox(current.first);//获取当前多边形边界框
            rtree.query(bgi::intersects(bbox), std::back_inserter(candidates));//通过rtree查找与当前边界框相交的多边形

            // 尝试合并所有相交的多边形
            for (const auto& candidate : candidates) {
                size_t j = candidate.second;//获取该多边形在inputs中的索引(第几个)
                if (merged[j] || i == j) continue;//跳过已经被处理过的环和自己本身

                if (polyIntersect(current.first, inputs[j].first)) {//精确检查是否相交
                    Polygon_t union_result = safe_union(current.first, inputs[j].first);
                    //执行多边形合并操作
                    // 检查合并后的多边形是否有效
                    if (bg::is_valid(union_result)) {//检查多边形是否有效
                        current.first = union_result;//更新当前多边形为合并后的结果
                        current.second.insert(current.second.end(),//？啥意思？没看懂，应该是合并索引，但是具体怎么合并的？
                            inputs[j].second.begin(),
                            inputs[j].second.end());
                        merged[j] = true;//标记为已处理
                    }
                }
            }

            combined_inner_polys.push_back(current);//将已经处理的环储存好(未处理的也会存)
        }

        // 第二步：检查合并后的多边形之间是否还有相交
        bool has_intersection = true;
        while (has_intersection) {
            has_intersection = false;//默认无相交
            std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> new_combined;
            std::vector<bool> processed(combined_inner_polys.size(), false);

            for (size_t i = 0; i < combined_inner_polys.size(); i++) {
                if (processed[i]) continue;

                auto current = combined_inner_polys[i];
                processed[i] = true;

                // 检查当前多边形是否与其他多边形相交
                for (size_t j = i + 1; j < combined_inner_polys.size(); j++) {
                    if (processed[j]) continue;

                    if (polyIntersect(current.first, combined_inner_polys[j].first)) {//检查
                        Polygon_t union_result = safe_union(current.first, combined_inner_polys[j].first);
                        //合并
                        if (bg::is_valid(union_result)) {
                            current.first = union_result;
                            current.second.insert(current.second.end(),
                                combined_inner_polys[j].second.begin(),
                                combined_inner_polys[j].second.end());
                            processed[j] = true;
                            has_intersection = true;
                        }
                    }
                }

                new_combined.push_back(current);
            }

            combined_inner_polys = std::move(new_combined);
        }

        // 第三步：最终检查，确保没有多边形相交
        for (size_t i = 0; i < combined_inner_polys.size(); i++) {
            for (size_t j = i + 1; j < combined_inner_polys.size(); j++) {
                if (polyIntersect(combined_inner_polys[i].first, combined_inner_polys[j].first)) {
                    // 尝试最后一次合并
                    Polygon_t union_result = safe_union(combined_inner_polys[i].first, combined_inner_polys[j].first);

                    if (bg::is_valid(union_result)) {
                        combined_inner_polys[i].first = union_result;
                        combined_inner_polys[i].second.insert(combined_inner_polys[i].second.end(),
                            combined_inner_polys[j].second.begin(),
                            combined_inner_polys[j].second.end());
                        combined_inner_polys.erase(combined_inner_polys.begin() + j);
                        j--; // 调整索引
                    }
                }
            }
        }

        // 创建合并后的内多边形并记录映射
        for (const auto& combined : combined_inner_polys) {
            if (combined.second.size() > 1) {
                Ring outer_ring = bg::exterior_ring(combined.first);
                if (bg::area(outer_ring) > area_threshold) {
                    RingNode node = formatNode(++maxRid, 1, outer_ring, orientation);
                    addNode(node);

                    std::vector<IdIndex> result_ids{ node.id };
                    MergeMap2 mm2{ IdIndex{0, 0}, combined.second, result_ids };
                    mergeMap2.push_back(mm2);

                    outputs.emplace_back(combined.first, node.id);
                }
            }
            else {
                outputs.emplace_back(combined.first, combined.second[0]);
            }
        }
    }


    std::vector<RingNode> FillBridge::compute_complex_polygon_merge(
        const std::vector<RingNode>& outers
        , const std::vector<RingNode>& inners
    ) {

        // 1. 准备内多边形
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> inner_polys;
        for (const auto& inner : inners) {
            Polygon_t poly;
            poly.outer() = inner.ring;//这是内多边形的边框
            bg::correct(poly);
            inner_polys.emplace_back(poly, std::vector<IdIndex>{inner.id});//其中包含了一个元素，就是这个多边形的键值
        }//精简参数，存档的就是多边形、键值。

        std::vector<std::pair<Polygon_t, IdIndex>> final_inners;
        polygonsMerge(inner_polys, final_inners, 1);//合并相交的内多边形，存在final——inners中

        //外多边形与内多边形合并
        std::vector<IdIndex> result_iis;//用来存储结果
        for (const auto& outer : outers) {//遍历外多边形
            Polygon_t outer_poly;
            outer_poly.outer() = outer.ring;
            bg::correct(outer_poly);//完成提取外多变形操作
            std::vector<Polygon_t> intersected_polys;//存放与当前外多边形相交的内多边形
            std::vector<IdIndex> intersected_iis;//以及其id

            // 与内多边形依次检查交集
            for (const auto& combined : final_inners) {//遍历刚才处理过的内多边形集合
                if (polyIntersect(outer_poly, combined.first)) {//检查是否相交
                    intersected_polys.emplace_back(combined.first);
                    intersected_iis.emplace_back(combined.second);//加入相交集合
                }
                else {
                    result_iis.emplace_back(combined.second);//不相交直接加入结果
                }
            }
            //去掉相交的内多边形
            final_inners.erase(
                std::remove_if(
                    final_inners.begin(),
                    final_inners.end(),
                    [&intersected_iis](const auto& pair) {
                        return std::find(
                            intersected_iis.begin(),
                            intersected_iis.end(),
                            pair.second
                        ) != intersected_iis.end();
                    }
                ),
                final_inners.end()
            );

            if (intersected_polys.size() > 0) {  // 有相交的内多边形
                MultiPolygon result;
                MultiPolygon intersected_mp;
                for (const auto& poly : intersected_polys) {//遍历有相交的内多边形
                    intersected_mp.push_back(poly);//存入临时内存
                }
                bg::difference(outer_poly, intersected_mp, result);//计算外多边形减去所有相交的内多边形的差集
                //类似于处理内外环互交形成子环的操作，子环存放在result中
                for (const auto& poly : result) {
                    Ring _ring = bg::exterior_ring(poly);//获取计算的差集的外环（还是以多边形进行操作）
                    if (bg::area(_ring) > area_threshold) {//检测新环面积是否过小，不是很小再进行下一步操作
                        RingNode node = formatNode(++maxRid, 1, _ring, -1);
                        addNode(node);
                        IdIndex _ii(node.id.id, 1);
                        std::vector<IdIndex> id_indices;
                        id_indices.emplace_back(_ii);
                        MergeMap2 mm2(outer.id, intersected_iis, id_indices);
                        mergeMap2.emplace_back(mm2);
                        result_iis.emplace_back(_ii);
                    }
                }
            }
            else {  // 与外多边形无相交
                result_iis.emplace_back(outer.id);
            }
        }

        //外多边形合并
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> outer_polys;
        for (IdIndex ii : result_iis) {
            RingNode rn = findNode(ii);
            Polygon_t poly;
            poly.outer() = rn.ring;
            bg::correct(poly);
            outer_polys.emplace_back(poly, std::vector<IdIndex>{ii});
        }
        std::vector<std::pair<Polygon_t, IdIndex>> final_polys;
        polygonsMerge(outer_polys, final_polys, -1);

        std::vector<RingNode> result;
        for (auto& pair : final_polys) {
            result.emplace_back(findNode(pair.second));
        }
        return result;
    }


    //生成环集
    void FillBridge::generateRings() {

        std::vector<RingNode> nodeHandles;//记录未处理的环集

        Ring ring = bg::exterior_ring(b_polygon);//取出外环并修正
        bg::correct(ring);
        RingNode node = formatNode(maxRid++, 1, ring, -1);//将外环包装为node对象格式，添加对应参数，maxrid是1
        nodeHandles.emplace_back(node);//将外环放到nodehandles中
        addNode(node);

        for (auto& inner : b_polygon.inners()) {//取出内环并修正
            bg::correct(inner);
            RingNode node = formatNode(maxRid++, 1, inner, 1);
            nodeHandles.emplace_back(node);
            addNode(node);//目的同外环
        }
        double t = 1;
        //遍历当前用于判断互交的节点集
        while (!nodeHandles.empty()) {
            std::vector<RingNode> outers;
            std::vector<RingNode> inners;
            for (const auto& rn : nodeHandles) {
                if (rn.orientation == -1) {//通过方向判断nodehandles中的单个环应该添加到outers中还是inners中
                    outers.emplace_back(rn);
                }
                else {
                    inners.emplace_back(rn);
                }
            }
            std::vector<RingNode> nodes = compute_complex_polygon_merge(
                outers, inners);
            std::vector<RingNode> offsetNodes;
            //偏移
            for (auto node : nodes) {
                std::vector<Ring> ring0 = offsetRing(node.ring, static_cast<double>(node.orientation) * offset * t, area_threshold);
                if (ring0.size() == 1) {
                    RingNode node0 = formatNode(node.id.id, ++node.id.index, ring0[0], node.orientation);
                    if (bg::area(node0.ring) > area_threshold) {
                        addNode(node0);
                        offsetNodes.emplace_back(node0);
                    }
                }
                else if (ring0.size() > 1) {
                    std::vector<IdIndex> id_indices;
                    for (auto r : ring0) {
                        RingNode node0 = formatNode(++maxRid, 1, r, node.orientation);
                        if (bg::area(node0.ring) > area_threshold) {
                            addNode(node0);
                            offsetNodes.emplace_back(node0);
                            IdIndex ii(maxRid, 1);
                            id_indices.emplace_back(ii);
                        }
                    }
                    //分裂映射
                    offsetMap.insert(std::pair<IdIndex, std::vector<IdIndex>>({ node.id.id,node.id.index }, id_indices));
                }
            }
            t = 1;
            //清空
            nodeHandles.clear();
            if (!offsetNodes.empty()) {
                if (offsetNodes.size() == 2) {
                    Polygon_t poly1, poly2;
                    bg::append(poly1.outer(), offsetNodes[0].ring);
                    bg::append(poly2.outer(), offsetNodes[1].ring);
                    std::vector<IdIndex> id_indices;
                    if (offsetNodes[0].orientation == 1 && isPolygonContained(poly1, poly2)) {
                        MergeMap mm(offsetNodes[0].id, offsetNodes[1].id, id_indices);
                        containMap.emplace_back(mm);
                        break;
                    }
                    if (offsetNodes[1].orientation == 1 && isPolygonContained(poly2, poly1)) {
                        MergeMap mm(offsetNodes[1].id, offsetNodes[0].id, id_indices);
                        containMap.emplace_back(mm);
                        break;
                    }
                }
                //重填
                nodeHandles.insert(nodeHandles.end(), offsetNodes.begin(), offsetNodes.end());
            }
        }
    }

    //形成树
    void FillBridge::formatTree() {
        //遍历不同环类型的节点集合
        for (auto& ringNode : ringNodes) {

            // 使用反向迭代器进行反序遍历
            std::vector<RingNode>& vec = ringNode.second;
            const auto it0 = ringNode.second.begin();
            if (it0->orientation == 1) { // 向外  倒序
                for (size_t i = ringNode.second.size() - 1; i > 0; i--) {
                    ringNode.second[i].children.emplace_back(ringNode.second[i - 1]);
                    ringNode.second[i - 1].parent = &ringNode.second[i];
                }
            }
            else { //向内 正序
                for (int i = 0; i < vec.size() - 1; i++) {
                    ringNode.second[i].children.emplace_back(ringNode.second[i + 1]);
                    ringNode.second[i + 1].parent = &ringNode.second[i];
                }
            }
        }
        //遍历分裂
        for (auto& map : offsetMap) { //偏移产生的分裂
            RingNode& it1 = findNode(map.first);
            for (const auto& ii : map.second) {
                RingNode& it2 = findNode(ii);
                if (it1.orientation == -1) {
                    it1.children.emplace_back(it2);
                    it2.parent = &it1;
                }
                else {
                    it1.parent = &it2;
                    it2.children.emplace_back(it1);
                }
            }
        }
        for (auto& map : containMap) {
            RingNode& it1 = findNode(map.ii1); //内
            RingNode& it2 = findNode(map.ii2);  //外
            it2.children.insert(it2.children.end(), it1.children.begin(), it1.children.end());
            it1.isHide = true;
        }
        for (auto& map : mergeMap2) {
            if (map.outer_ii.id == 0 && map.outer_ii.index == 0) {  //内多边形合并
                for (const auto& merged : map.merged_iis) {
                    RingNode& merge_rn = findNode(merged);
                    for (const auto& inner_ii : map.inner_iis) {
                        RingNode& inner_rn = findNode(inner_ii);
                        merge_rn.children.insert(merge_rn.children.end(),
                            inner_rn.children.begin(), inner_rn.children.end());
                    }
                }
            }
            else {  //外多边形  合并
                RingNode& outer = findNode(map.outer_ii);
                for (const auto& merged : map.merged_iis) {
                    RingNode& merge_rn = findNode(merged);
                    size_t size = ringNodes.size();

                    outer.parent->children.emplace_back(merge_rn);
                    merge_rn.parent = outer.parent;
                    IdIndex _ii = outer.id;
                    //将outer从父节点的子节点集中移除
                    outer.parent->children.erase(std::remove_if(outer.parent->children.begin(), outer.parent->children.end(),
                        [_ii](RingNode x) { return x.id.id == _ii.id && x.id.index == _ii.index; }),
                        outer.parent->children.end());
                    //将内多边形的子节点作为merge_rn的子节点
                    for (const auto& inner_ii : map.inner_iis) {
                        RingNode& inner_rn = findNode(inner_ii);
                        merge_rn.children.insert(merge_rn.children.end(),
                            inner_rn.children.begin(), inner_rn.children.end());
                    }
                }

            }
        }

    }
    //树遍历
    void FillBridge::dfs(
        RingNode& node,
        std::vector<IdIndex>& visited
    ) {
        for (auto& one : visited) {
            if (one == node.id) {
                return;
            }
        }
        visited.emplace_back(node.id);
        for (auto& child : node.children) {
            dfs(findNode(child.id), visited);
            handleBridge(node.id, child.id);
        }
    }

    // 计算点到线段的最近点
    Point_t FillBridge::closest_point_on_segment(const Point_t& p, const Point_t& seg_start, const Point_t& seg_end) {
        const double x = bg::get<0>(p);
        const double y = bg::get<1>(p);
        const double x1 = bg::get<0>(seg_start);
        const double y1 = bg::get<1>(seg_start);
        const double x2 = bg::get<0>(seg_end);
        const double y2 = bg::get<1>(seg_end);

        // 线段向量
        const double dx = x2 - x1;
        const double dy = y2 - y1;

        // 如果线段长度为零，返回起点
        if (dx == 0 && dy == 0) {
            return seg_start;
        }

        // 计算投影参数t
        double t = ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy);

        // 限制 t 在 [0,1] 范围内，确保投影点在线段上
        t = std::max<double>(0.0, std::min<double>(1.0, t));

        // 计算投影点坐标
        return Point_t(bg::get<0>(seg_start) + t * dx,
            bg::get<1>(seg_start) + t * dy);
    }

    // 在环的边上查找离给定点最近的点
    Point_t FillBridge::find_closest_point_on_ring_edges(Ring& ring, const Point_t& p0, size_t& e_index1) {
        Point_t closest = ring[0]; // 初始化为Ring的第一个点
        double minDist = bg::distance(p0, closest);

        // 遍历Ring的每条边
        for (size_t i = 0; i < ring.size(); ++i) {
            // 获取当前边的两个端点
            Point_t a = ring[i];
            Point_t b = ring[(i + 1) % ring.size()]; // 处理闭合边

            // 计算点到边的最近点
            Point_t projection = closest_point_on_segment(p0, a, b);

            // 计算距离
            double dist = bg::distance(p0, projection);

            // 更新最近点
            if (dist < minDist) {
                minDist = dist;
                closest = projection;
                e_index1 = i;
            }
        }
        return closest;
    }

    Point_t FillBridge::find_point_at_distance_clockwise(Ring& ring, const Point_t& start_point,
        size_t _index, double d, size_t& e_index)
    {
        // 确保环是闭合的（首尾点相同）
        Ring normalized_r1 = ring;
        double remaining_distance = d;
        size_t current_edge = _index;
        Point_t current_point = start_point;
        while (remaining_distance > 0) {
            // 获取当前边的终点
            Point_t next_point = normalized_r1[(current_edge + 1) % normalized_r1.size()];

            // 计算当前边剩余长度
            double edge_length = bg::distance(current_point, next_point);

            if (edge_length >= remaining_distance) {
                // 目标点在当前边上
                double ratio = remaining_distance / edge_length;
                double x = current_point.x() + ratio * (next_point.x() - current_point.x());
                double y = current_point.y() + ratio * (next_point.y() - current_point.y());
                e_index = current_edge;
                return Point_t(x, y);
            }
            else {
                // 目标点在下一条边上，更新剩余距离
                remaining_distance -= edge_length;
                current_point = next_point;
                current_edge = (current_edge + 1) % normalized_r1.size();
            }
        }
        e_index = _index;
        return start_point;
    }


    //点在环上的位置索引
    int FillBridge::findPointIndex(const Ring& ring, const Point_t& p0) {
        for (size_t i = 0; i < ring.size(); ++i) {
            if (equal(ring[i], p0)) {
                return static_cast<int>(i);
            }
        }
        return -1; // 未找到
    }

    bool FillBridge::does_segment_cross_ring(const Segment& seg, const Ring& ring) {
        // 1. 获取线段端点
        const Point_t& p0 = seg.first;
        const Point_t& p1 = seg.second;

        // 2. 遍历环的所有边
        for (size_t i = 0; i < ring.size(); i++) {
            // 当前边的两个端点
            const Point_t& r0 = ring[i];
            const Point_t& r1 = ring[(i + 1) % ring.size()];

            // 3. 排除线段端点与环顶点重合的情况
            if (equal(p0, r0) || equal(p0, r1) ||
                equal(p1, r0) || equal(p1, r1)) {
                continue;
            }

            // 4. 检查线段与边是否相交
            if (bg::intersects(seg, Segment(r0, r1))) {
                std::vector<Point_t> intersection_points;
                bg::intersection(seg, Segment(r0, r1), intersection_points);
                if (intersection_points.size() > 1) {
                    return true;
                }
            }
        }
        return false;
    }


    void FillBridge::handleBridge(IdIndex o_ii, IdIndex i_ii)
    {   //判断环上是否存在桥接
        for (auto& b : bridges) {
            if (b.to_ii == i_ii) {
                return;
            }
        }
        RingNode& inner = findNode(i_ii);
        RingNode& outer = findNode(o_ii);

        double _offset = offset;
        Point_t p0, p1, p2, p3;
        size_t p0_index = 0, p1_index = 0, p2_index = 0, p3_index = 0, _index = 0;
        double perimeter = 0;
        //内环所有边的中点（降序）及周长
        std::vector<MidPoints> points = find_mid_points_on_ring(inner.ring, perimeter);
        if (bridges.size() == 0) {
            p0 = points[0].mid_point;
            p0_index = points[0].index;
        }
        else {
            BridgeMap bridge = bridges.back();
            size_t _index = findPointIndex(inner.ring, bridge.from);
            find_further_point(points, bridge.from, _index, p0, p0_index);
        }
        
        bool findBridge = false;
        double total_len = 0;
        int idx = 0;
        const double min_distance = 0.95 * _offset; // 最小距离阈值
        const int max_attempts = 200; // 最大尝试次数限制
        int attempt_count = 0;

        while (!findBridge && attempt_count < max_attempts) {
            attempt_count++;

            p1 = find_point_at_distance_clockwise(inner.ring, p0, p0_index, _offset, p1_index);
            //得到离内环p0点最近的外环点p2
            p2 = find_closest_point_on_ring_edges(outer.ring, p0, p2_index);
            p3 = find_closest_point_on_ring_edges(outer.ring, p1, p3_index);
            total_len += _offset;

            // 检查点对距离是否过近
            if (bg::distance(p0, p1) < min_distance 
                || bg::distance(p2, p3) < min_distance
                || bg::distance(p3, p2) < min_distance
                ) {
                if (total_len > perimeter) {
                    //_offset *= 0.9;
                    total_len = 0;
                    idx = (idx + 1) % points.size();
                    p0 = points[idx].mid_point;
                    p0_index = points[idx].index;
                }
                else {
                    p0 = p1;
                    p0_index = p1_index;
                }
                continue;
            }

            Segment seg1(p2, p1), seg2(p3, p0), seg3(p2, p0), seg4(p3, p1);
            std::vector<Point_t> temp, temp1;
            bg::intersection(seg1, seg2, temp);
            bg::intersection(seg3, seg4, temp1);
            if (temp.empty() || temp.size() > 1) {
                // 处理无交点情况
                if (total_len > perimeter) {
                    _offset *= 0.9;
                    total_len = 0;
                    idx = (idx + 1) % points.size();
                    p0 = points[idx].mid_point;
                    p0_index = points[idx].index;
                }
                else {
                    p0 = p1;
                    p0_index = p1_index;
                }
                continue;
            }

             // 验证桥接线段长度
            if (bg::distance(p1, p2) > 1.5 * _offset || bg::distance(p0, p3) > 1.5 * _offset) {
                // 处理过长线段
                if (total_len > perimeter) {
                    _offset *= 0.9;
                    total_len = 0;
                    idx = (idx + 1) % points.size();
                    p0 = points[idx].mid_point;
                    p0_index = points[idx].index;
                }
                else {
                    p0 = p1;
                    p0_index = p1_index;
                }
                continue;
            }
            // 检查新桥接是否与已有桥接交叉
            bool crossesExisting = false;
            for (const auto& existing_bridge : bridges) { 
                // 检查新桥接线段是否与已有桥接线段相交
                if (bg::distance(p0, existing_bridge.from) < _offset 
                    || bg::distance(p1, existing_bridge.from2) < _offset
                    || bg::distance(p0, existing_bridge.from2) < _offset
                    || bg::distance(p1, existing_bridge.from) < _offset
                    ) {
                    crossesExisting = true;
                    break;
                }
            }

            if (crossesExisting) {
                // 如果与已有桥接交叉，尝试调整位置
                if (total_len > perimeter) {
                    _offset *= 0.9;
                    total_len = 0;
                    idx = (idx + 1) % points.size();
                    p0 = points[idx].mid_point;
                    p0_index = points[idx].index;
                }
                else {
                    p0 = p1;
                    p0_index = p1_index;
                }
                continue;
            }
            //// 插入点并记录桥接
            if (p0_index < p1_index) {
                insert_point_on_ring(inner.ring, p0, p0_index);
                insert_point_on_ring(inner.ring, p1, p1_index);
            }
            else {
                insert_point_on_ring(inner.ring, p1, p1_index);
                insert_point_on_ring(inner.ring, p0, p0_index);
            }
            if (p2_index < p3_index) {
                insert_point_on_ring(outer.ring, p2, p2_index);
                insert_point_on_ring(outer.ring, p3, p3_index);
            }
            else {
                insert_point_on_ring(outer.ring, p3, p3_index);
                insert_point_on_ring(outer.ring, p2, p2_index);
            }
            /*insertPointIntoRing(inner.ring, p0);
            insertPointIntoRing(inner.ring, p1);
            insertPointIntoRing(outer.ring, p2);
            insertPointIntoRing(outer.ring, p3);*/

            bridges.emplace_back(p2, p1, p3, p0, o_ii, i_ii);
            findBridge = true;
        }
        if (attempt_count >= max_attempts) {
            // 处理无法找到合适桥接的情况
            std::cout << "Warning: Could not find non-crossing bridge after "
                << max_attempts << " attempts." << std::endl;
        }
    }


    int FillBridge::findIndex(std::vector<Point_t> points, const Point_t& p0) {
        for (int i = 0; i < points.size(); ++i) {
            if (equal(points[i], p0)) {
                return static_cast<int>(i);
            }
        }
        return -1; // 未找到
    }


     //递归遍历环
    void FillBridge::traverseRing(
        RingNode& node,
        Point_t& start,
        Point_t& end,
        bool isOutermostLayer) {

        int size = node.ring.size();
        int s_index = findPointIndex(node.ring, start);
        int e_index = findPointIndex(node.ring, end);
        std::vector<Point_t> c_points;
        std::vector<Point_t> cc_points;
        std::vector<Point_t> points;
        for (int i = s_index; !equal(node.ring[i], end); i = (i + 1) % size) {
            c_points.emplace_back(node.ring[i]);
        }
        c_points.emplace_back(end);
        for (int i = s_index; !equal(node.ring[i], end); i = (i - 1 + size) % size) {
            cc_points.emplace_back(node.ring[i]);
        }
        cc_points.emplace_back(end);

        if (cc_points.size() > c_points.size()) {
            points.assign(cc_points.begin(), cc_points.end());
        }
        else
        {
            points.assign(c_points.begin(), c_points.end());
        }
        int index = 0;
        do {
            path.emplace_back(points[index]);
            BridgeMap* bridge = nullptr;
            //判断该点是否为桥接点
            for (auto& bm : bridges) {
                if (equal(bm.from, points[index]) || equal(bm.from2, points[index])) {
                    bridge = &bm;
                    break;
                }
            }
            if (bridge != nullptr) {
                RingNode& rn = findNode(bridge->to_ii); //to 环
                if (equal(bridge->from, points[index])) {
                    traverseRing(rn, bridge->to, bridge->to2, false);
                    index = findIndex(points, bridge->from2);
                }
                else if (equal(bridge->from2, points[index])) {
                    traverseRing(rn, bridge->to2, bridge->to, false);
                    index = findIndex(points, bridge->from);
                }
                if (index == -1) {
                    std::cout << bridge->from.x() << "," << bridge->from.y() << std::endl;
                    break;
                }
                path.emplace_back(points[index]);
            }
            index += isOutermostLayer ? -1 : 1;
            // 额外检查：防止index溢出导致死循环
            if (index > points.size() || index < 0) {
                break;
            }
        } while (index > 0 && index < points.size());
    }
} // namespace Slic3r


