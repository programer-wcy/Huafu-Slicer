//
// Created by ylzha on 2025/2/25.
//
#ifndef slic3r_FillBridge_hpp_
#define slic3r_FillBridge_hpp_

#include "../libslic3r.h"
#include <boost/geometry.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include "FillBase.hpp"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Slic3r{

class Surface;

// 定义点类型
typedef bg::model::d2::point_xy<double> Point_t;

//定义线段
typedef bg::model::segment<Point_t> Segment;

// 定义多边形类型
typedef bg::model::polygon<Point_t> Polygon_t;

// 多边形集合
typedef bg::model::multi_polygon<Polygon_t> MultiPolygon;

// 定义环类型
typedef bg::model::ring<Point_t> Ring;

// 定义折线类型（用于表示直线）
typedef bg::model::linestring<Point_t> Linestring;
typedef bg::model::box<Point_t> Box_t;

// 2. 图结构定义（Boost.Graph）
typedef boost::adjacency_list<
    boost::vecS,               // 边存储：动态数组
    boost::vecS,               // 顶点存储：动态数组
    boost::undirectedS,        // 无向图（桥接点连接无方向）
    Point_t,                   // 顶点属性：存储 Point_t（坐标）
    boost::no_property         // 边属性：无额外信息
> ShapeGraph;

typedef boost::graph_traits<ShapeGraph>::edge_iterator EdgeIter;  // 边迭代器
typedef boost::graph_traits<ShapeGraph>::vertex_iterator VertexIter; // 顶点迭代器
typedef boost::graph_traits<ShapeGraph>::vertex_descriptor Vertex; // 顶点句柄
typedef boost::graph_traits<ShapeGraph>::edge_descriptor   Edge;   // 边句柄

// 定义路径类型（顶点序列）
typedef std::vector<Vertex> Path;

struct PointCompare {
    bool operator()(const Point_t& a, const Point_t& b) const {
        if (a.x() != b.x()) return a.x() < b.x();
        return a.y() < b.y();
    }
};

struct MidPoints {
    size_t index; //边索引
    Point_t mid_point; //中点
    double distance; //某点离边的距离
    MidPoints(size_t index, Point_t mid_point, double distance) :index(index)
        , mid_point(mid_point), distance(distance) {
    }
};

struct IdIndex {
    size_t id{};
    size_t index{};
    // 定义小于运算符
    bool operator<(const IdIndex& other) const {
        return id < other.id || (id == other.id && index < other.index);
    }
    IdIndex() = default;
    IdIndex(const size_t id, const size_t index) : id(id), index(index) {}

    bool operator==(const IdIndex& ii) const {
        return id == ii.id && index == ii.index;
    }
};

//环路径树结构中节点类型
struct RingNode {
    IdIndex id;
    Ring ring; //当前环
    int orientation; //方向 -1 向内 1 向外
    std::vector<RingNode> children; //子节点集
    RingNode* parent; //父节点
    bool isHide{ false }; //在最终路径中是否呈现 false 呈现  true 不呈现
    // 构造函数，方便初始化
    RingNode(const IdIndex id, const Ring& ring, const int orientation = 0)
        : id(id), ring(ring), orientation(orientation), parent(nullptr) {
    }

    // 重载==运算符方便比较
    bool operator==(const RingNode& other) const {
        return id.id == other.id.id && id.index == other.id.index;
    }

};

//桥接映射
struct BridgeMap {
    Point_t from;
    Point_t to;
    Point_t from2;
    Point_t to2;
    IdIndex from_ii;
    IdIndex to_ii;
    // 默认构造函数
    BridgeMap()
        : from(0, 0), to(0, 0), from2(0, 0), to2(0, 0),  // 给 Point_t 传默认参数
        from_ii({}), to_ii({}) {
    }                       

    BridgeMap(Point_t from, Point_t to, 
        Point_t from2, Point_t to2, IdIndex from_ii
        , IdIndex to_ii
    ) :
        from(from), to(to), 
        from2(from2), to2(to2), 
        from_ii(from_ii), to_ii(to_ii){
    }
};

//合并映射
struct MergeMap {
    IdIndex ii1; //环1
    IdIndex ii2; //环2
    std::vector<IdIndex> nodes; //产生的环集
    // 构造函数
    MergeMap(const IdIndex& i1, const IdIndex& i2, const std::vector<IdIndex>& nodeList)
        : ii1(i1), ii2(i2), nodes(nodeList) {
    }
};

struct MergeMap2 {
    IdIndex outer_ii; //外多边形
    std::vector<IdIndex> inner_iis; //多个内多边形
    std::vector<IdIndex> merged_iis; //合并产生的多边形
    MergeMap2(const IdIndex& outer_ii, const std::vector<IdIndex>& inner_iis,
        const std::vector<IdIndex> merged_iis)
        :outer_ii(outer_ii), inner_iis(inner_iis), merged_iis(merged_iis) {
    }
};

struct RingVertices {
    IdIndex ii;  //环ID
    std::vector<std::pair<Point_t, Point_t>> bridge_pairs;  //桥接点对
    std::vector<Vertex> ring_vertices; //环顶点集
    RingVertices(const IdIndex& ii, const std::vector<std::pair<Point_t, Point_t>>& bridge_pairs,
        const std::vector<Vertex>& ring_vertices) :
        ii(ii), bridge_pairs(bridge_pairs), ring_vertices(ring_vertices) {
    }
};

class FillBridge : public Fill {

public:
    ~FillBridge() override = default;
    bool is_self_crossing() override { return false; }
//protected:
    Fill* clone() const override { return new FillBridge(*this); }
    Polylines fill_surface(const Surface* surface, const FillParams& params) override;

private:
    std::map<size_t, std::vector<RingNode>> ringNodes;
    std::vector<BridgeMap> bridges;
    std::map<IdIndex, std::vector<IdIndex>> offsetMap;
    std::vector<MergeMap> containMap;
    std::vector<MergeMap2> mergeMap2;
    size_t maxRid = 1;
    std::vector<Point_t> path;
    bool isFront = false;
    double offset = 80000; //偏移距离
    
    double area_threshold = sqr(scaled<double>(0.5)); //面积阈值
    Polygon_t b_polygon;
    Polylines b_polylines;

    ShapeGraph graph;
    std::vector<Vertex> vertices; // 存储所有顶点句柄
    std::vector<RingVertices> rv_vertices; //存储环上顶点

private:

    // 判断两个环是否相交（不包括包含关系）
    bool polyIntersect(const Polygon_t poly1, const Polygon_t poly2);

    //获取环偏移后的环（可能为空）
    std::vector<Ring> offsetRing(const Ring& ring, double distance, double area_threshold);

    //生成环集
    void generateRings();
    //形成节点
    RingNode formatNode(size_t id, size_t index, const Ring& ring, int orientation);
    void addNode(RingNode& node);
    void removeNode(IdIndex id);
    //查找节点
    RingNode& findNode(IdIndex ii);

    //形成树
    void formatTree();
    // 深度优先搜索遍历树
    void dfs(RingNode& node,std::vector<IdIndex>& visited);
    Polygon_t safe_union(const Polygon_t& a, const Polygon_t& b);

    // 在环上按顺时针方向查找距离给定点d的另一个点
    Point_t find_point_at_distance_clockwise(Ring& ring, const Point_t& start_point, 
        size_t _index, double d, size_t& e_index);

    bool equal(Point_t p1, Point_t p2);

 
     void handleBridge(IdIndex o_ii, IdIndex i_ii);

    //递归遍历环
    //void traverseRing(BridgeMap bm, bool isOutermostLayer );
    void traverseRing(
        RingNode& node,
        Point_t& start,
        Point_t& end,
        bool isOutermostLayer);
    std::vector<RingNode> compute_complex_polygon_merge(
        const std::vector<RingNode>& outers
        , const std::vector<RingNode>& inners
    );
    bool isPolygonContained(const Polygon_t& polyA, const Polygon_t& polyB);

    Point_t closest_point_on_segment(const Point_t& p, const Point_t& seg_start, const Point_t& seg_end);
    Point_t find_closest_point_on_ring_edges(Ring& ring, const Point_t& p0, size_t& e_index1);
    int findPointIndex(const Ring& ring, const Point_t& p0);
    bool does_segment_cross_ring(const Segment& seg, const Ring& ring);

    int findIndex(std::vector<Point_t> points, const Point_t& p0);
    auto get_bbox(const Polygon_t& poly);
    void polygonsMerge(
        std::vector<std::pair<Polygon_t, std::vector<IdIndex>>> inputs,
        std::vector<std::pair<Polygon_t, IdIndex>>& outputs, size_t orientation);


    //环所有边按照边长遍历边的中点
    std::vector<MidPoints> find_mid_points_on_ring(const Ring& ring, double& perimeter) {
        // 存储边长度和索引的结构
        std::vector<MidPoints> points;
        perimeter = 0;
        // 遍历环的边（注意闭合环最后一点与起点相连）
        for (size_t i = 0; i < ring.size(); ++i) {
            const Point_t& p1 = ring[i];
            const Point_t& p2 = ring[(i + 1) % ring.size()]; // 处理闭合边
            // 计算线段中点
            Point_t midpoint(
                (bg::get<0>(p1) + bg::get<0>(p2)) / 2,
                (bg::get<1>(p1) + bg::get<1>(p2)) / 2
            );
            double len = bg::distance(p1, p2);
            MidPoints m(i, midpoint, len);
            points.push_back(m);
            perimeter += len;
        }
        //根据距离降序排序
        std::sort(points.begin(), points.end(),
            [](const MidPoints& a, const MidPoints& b) {
                return a.distance > b.distance; // 降序：a的distance大于b时，a排在前面
            });
        return points;

    }

    //距离点point最远的长边中点
    void find_further_point(std::vector<MidPoints> midpoints, 
        Point_t point, size_t _index,Point_t& point2, size_t& point2_index) {
        double max_distance = 0.0;
        for (auto& mp : midpoints) {
            double d = bg::distance(mp.mid_point, point);
            if (max_distance < d && mp.index != _index) {
                point2 = mp.mid_point;
                point2_index = mp.index;
                max_distance = d;
            }
        }
    }



    void insertPointIntoRing(Ring& r, const Point_t& p0) {
        if (r.empty()) {
            r.push_back(p0);
            return;
        }

        if (r.size() == 1) {
            r.push_back(p0);
            return;
        }

        std::size_t n = r.size();
        double min_distance = (std::numeric_limits<double>::max)();
        std::size_t insert_index = 0;

        for (std::size_t i = 0; i < n; ++i) {
            const Point_t& p1 = r[i];
            const Point_t& p2 = r[(i + 1) % n];

            // 创建线段
            Segment seg(p1, p2);

            // 计算点到线段的距离
            double dist = bg::distance(p0, seg);

            if (dist < min_distance) {
                min_distance = dist;
                insert_index = i + 1;
            }
        }

        // 在找到的位置插入点
        if (insert_index == n) {
            r.push_back(p0);
        }
        else {
            r.insert(r.begin() + insert_index, p0);
        }
    }


    /**
     * 根据Point_t查找对应的Vertex（考虑浮点精度）
     * @param graph 图结构
     * @param p0 目标点
     * @param epsilon 精度阈值（默认1e-6，处理浮点数误差）
     * @return 所有匹配的Vertex（可能多个，若存在重复点）
     */
    std::vector<Vertex> findVertexByPoint(const ShapeGraph& graph, const Point_t& p0, double epsilon = 1e-6) {
        std::vector<Vertex> matchedVertices;

        // 获取所有顶点的迭代器范围
        std::pair<VertexIter, VertexIter> vertices = boost::vertices(graph);

        // 遍历所有顶点，比较坐标
        for (VertexIter it = vertices.first; it != vertices.second; ++it) {
            Vertex v = *it;
            const Point_t& p = graph[v]; // 获取当前顶点的坐标属性
            if (equal(p, p0)) {
                matchedVertices.push_back(v);
            }

            //// 计算两点距离（Boost.Geometry提供的距离函数）
            //double distance = bg::distance(p, p0);

            //// 若距离小于精度阈值，视为匹配
            //if (distance < epsilon) {
            //    matchedVertices.push_back(v);
            //}
        }

        return matchedVertices;
    }


    // 6. 用Hierholzer算法生成欧拉路径
    std::vector<Vertex> findEulerPath(ShapeGraph graph, Vertex start) {
        std::vector<Vertex> path;
        std::stack<Vertex> currentPath;
        currentPath.push(start);

        while (!currentPath.empty()) {
            Vertex u = currentPath.top();

            // 检查当前顶点是否有未访问的边
            auto [eBegin, eEnd] = boost::out_edges(u, graph);
            if (eBegin != eEnd) {
                // 取一条未访问的边，移动到相邻顶点
                Edge e = *eBegin;
                Vertex v = boost::target(e, graph);
                boost::remove_edge(e, graph); // 移除边（标记为已访问）
                currentPath.push(v);
            }
            else {
                // 无未访问边，加入路径并回溯
                path.push_back(u);
                currentPath.pop();
            }
        }

        std::reverse(path.begin(), path.end()); // 反转得到正确的路径顺序
        return path;
    }


    void insert_point_on_ring(Ring& ring, Point_t& point, size_t index) {
        // 创建环的副本进行操作，而不是修改原始数据
        Ring ring_copy = ring;
        const auto insertPos = ring_copy.begin() + index + 1;
        ring_copy.insert(insertPos, point);
        bg::correct(ring_copy); //调整方向为顺时针
        // 更新原始环
        ring = std::move(ring_copy);
    }

     
    /**
     * 辅助函数：将点转换为可哈希的字符串（用于检测非连续重复）
     * 处理浮点数精度问题（保留6位小数）
     */
    std::string point_to_string(const Point_t& p) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6);  // 控制精度，避免浮点数误差导致的误判
        ss << bg::get<0>(p) << "," << bg::get<1>(p);
        return ss.str();
    }

    /**
     * 改进版：去除环中所有重复点（包括非连续重复）
     * 保留环的闭合性（首尾点相同）
     */
    void remove_duplicate_points(Ring& input_ring) {
        if (input_ring.size() <= 1) {
            return;  // 空环或单点环无需处理
        }

        Ring result;
        std::unordered_set<std::string> seen_points;  // 记录已添加的点（字符串哈希）

        // 处理第一个点（环的起点）
        const Point_t& first_point = input_ring[0];
        std::string first_str = point_to_string(first_point);
        result.push_back(first_point);
        seen_points.insert(first_str);

        // 遍历中间点（跳过第一个和最后一个，最后一个单独处理以保持闭合性）
        size_t last_idx = input_ring.size() - 1;
        for (size_t i = 1; i < last_idx; ++i) {
            const Point_t& current = input_ring[i];
            std::string current_str = point_to_string(current);

            // 如果当前点未出现过，则添加
            if (seen_points.find(current_str) == seen_points.end()) {
                result.push_back(current);
                seen_points.insert(current_str);
            }
        }

        // 处理最后一个点（确保环闭合）
        const Point_t& last_point = input_ring[last_idx];
        std::string last_str = point_to_string(last_point);

        // 环的最后一个点应该与第一个点相同（闭合性），无论是否重复都保留
        // 若原环最后一个点与第一个点不同，则检查是否为新点
        if (bg::equals(last_point, first_point)) {
            result.push_back(last_point);  // 保持闭合
        }
        else {
            if (seen_points.find(last_str) == seen_points.end()) {
                result.push_back(last_point);
            }
            // 若最后一个点是重复点且不与起点相同，则不添加（避免破坏闭合性）
        }

        // 特殊情况：如果处理后只剩一个点，手动闭合
        if (result.size() == 1) {
            result.push_back(result[0]);
        }

        input_ring = std::move(result);
    }

    // 查找函数
    Vertex find_vertex_by_point(const std::vector<Vertex>& ring_vertices,
        const ShapeGraph& graph,
        const Point_t& p) {
        for (Vertex v : ring_vertices) {
            if (bg::equals(graph[v], p)) {
                return v;
            }
        }
        return boost::graph_traits<ShapeGraph>::null_vertex();
    }

 

    // 判断i是否在_vertex_list中
    bool isInVertexList(int i, const std::vector<int>& _vertex_list) {
        // 使用std::find查找元素i
        auto it = std::find(_vertex_list.begin(), _vertex_list.end(), i);

        // 如果迭代器不等于end()，说明找到元素
        return it != _vertex_list.end();
    }

    /**
     * 删除图中所有度数为1的顶点及其关联的边
     * @param graph 要操作的图
     */
    void removeDegreeOneVertices(ShapeGraph& graph) {
        // 1. 收集所有度数为1的顶点
        std::vector<Vertex> verticesToRemove;

        // 获取顶点迭代器范围
        std::pair<VertexIter, VertexIter> vertices = boost::vertices(graph);

        // 遍历所有顶点，筛选出度数为1的顶点
        for (VertexIter it = vertices.first; it != vertices.second; ++it) {
            Vertex v = *it;
            // 计算顶点的度数（连接的边数）
            std::size_t degree = boost::degree(v, graph);

            if (degree == 1) {
                verticesToRemove.push_back(v);
            }
        }

        // 2. 删除所有收集到的顶点（会自动删除关联的边）
        for (Vertex v : verticesToRemove) {
            // 注意：使用vecS作为顶点存储时，删除顶点后可能影响其他顶点的描述符
            // 但这里我们使用提前收集的顶点列表，避免迭代器失效问题
            boost::remove_vertex(v, graph);
        }
    }

};
}; // namespace Slic3r
#endif // slic3r_FillBridge_hpp_
