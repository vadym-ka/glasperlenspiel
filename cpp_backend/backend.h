#include <iostream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <random>
#include <assert.h>
#include <math.h>
#include <memory>

#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"


namespace {

using Id = size_t;

template<class T1, class T2>
std::ostream& operator<< (std::ostream &out, std::pair<T1, T2> pair) { return out << "(" << pair.first << ", " << pair.second << ")";}
template<class T> std::ostream& operator<<(std::ostream& out, std::vector<T> vec) { out<<"("; for (auto& v: vec) out<<v<<", "; return out<<")"; }
template<class T> std::ostream& operator<<(std::ostream& out, std::set<T> vec) { out<<"("; for (auto& v: vec) out<<v<<", "; return out<<")"; }
template<class L, class R> std::ostream& operator<<(std::ostream& out, std::map<L, R> vec) { out<<"("; for (auto& v: vec) out<<v<<", "; return out<<")"; }

std::mt19937 rnd(19937 + 228 + 14 / 88 + 239 + 24111997 % 322);

}  // namespace


namespace euclid_plane {

typedef long double ld;

static constexpr ld kEps = 1e-9;

namespace {

ld hypotenuse_length(ld x, ld y) {
  x = abs(x);
  y = abs(y);
  if (x > y) std::swap(x, y);
  if (abs(y) < kEps) return 0;
  ld d = x / y;
  return y * sqrt(1 + d * d);
}

ld determinant(ld a, ld b, ld c, ld d) {
  return a * d - b * c;
}

}  // namespace

// Point, Line, Circle
// Relations between objects
// Tools to build graph hmmmmm, toolset : [LineByTwoPoint, PointByTwoLines, ]
// Graph view :) (allow us not to delete objects from graph)
// Visualisation with Geogebra (just generate js code and save pictures!!!)
// What is observation space? (hm, jsoned graph!!!!!!!)

// Tools :
// args,
// requirements (constraints)
// result (unique?)

enum class Type {
  Point, Line, Circle, None
};

struct Object {
  virtual Type type() const { return Type::None; }
  virtual bool operator==(const Object& other) const { assert(false); }
};

struct Point : Object {
  Point(ld x, ld y) : x(x), y(y) {}
  Point() : Point(0, 0) {}
  Type type() const override { return Type::Point; }
  bool operator<(const Point& other) const {
    if (abs(x - other.x) > kEps) return x < other.x;
    if (abs(y - other.y) > kEps) return y < other.y;
    return false;
  }
  bool operator==(const Point& other) const {
    if (abs(x - other.x) > kEps) return false;
    if (abs(y - other.y) > kEps) return false;
    return true;
  }
  bool operator!=(const Point& other) const { return !(*this == other); }
  // bool operator==(const Point& other) const { return abs(x - other.x) < kEps && abs(y - other.y) < kEps; }
  Point& operator+=(const Point& other) { x += other.x; y += other.y; return *this; }
  Point operator+(const Point& other) const { Point ans = *this; return ans += other; }
  Point& operator-=(const Point& other) { x -= other.x; y -= other.y; return *this; }
  Point operator-(const Point& other) const { Point ans = *this; return ans -= other; }
  ld operator*(const Point& other) const { return x * other.x + y * other.y; }
  ld operator^(const Point& other) const { return x * other.y - y * other.x; }


  ld x, y;
};
std::ostream& operator<<(std::ostream& out, const Point& p) { return out << "( " << p.x << ", " << p.y << " )"; }

struct Line : Object {
  Line(ld a, ld b, ld c) : a(a), b(b), c(c) {}

  Line() : Line(0, 0, 0) {}
  Type type() const override { return Type::Line; }

  bool operator==(const Line &other) const {
    Line l1 = this->normed();
    Line l2 = other.normed();
    return abs(l1.a - l2.a) < kEps &&
           abs(l1.b - l2.b) < kEps &&
           abs(l1.c - l2.c) < kEps;
  }

  bool operator!=(const Line &other) const { return !(*this == other); }

  Line normed() const {
    Line l = *this;
    ld norm_coef = hypotenuse_length(a, b);
    l.a /= norm_coef;
    l.b /= norm_coef;
    l.c /= norm_coef;
    if (l.a < -kEps || (abs(l.a) < kEps && l.b < -kEps)) {
      l.a *= -1;
      l.b *= -1;
      l.c *= -1;
    }
    return l;
  }

  bool contains(Point p) const { return abs(a * p.x + b * p.y + c) < kEps; }

  ld a, b, c;
};

std::ostream& operator<<(std::ostream& out, const Line& l) { return out << "( " << l.a << ", " << l.b << ", " << l.c << " )"; }

struct Circle : Object {
  Circle(const Point &p, ld radius);
  Point center;
  ld radius;
};

namespace tools {

Line line_by_two_points(const Point &p1, const Point &p2) {
  assert(p1 != p2);
  ld a = p1.y - p2.y;
  ld b = p2.x - p1.x;
  ld c = -a * p1.x - b * p1.y;
  Line l = Line(a, b, c);
  return l.normed();
}

bool are_parallel(const Line &l1, const Line &l2) {
  ld det = determinant(l1.a, l1.b, l2.a, l2.b);
  return abs(det) < kEps;
}

bool are_perpendicular(const Line &l1, const Line &l2) {
  ld dot = l1.a * l2.a + l1.b * l2.b;
  return abs(dot) < kEps;
}

Point intersect_two_lines(const Line &l1, const Line &l2) {
  assert(!are_parallel(l1, l2));
  ld det = determinant(l1.a, l1.b, l2.a, l2.b);
  Point p;
  p.x = -determinant(l1.c, l1.b, l2.c, l2.b) / det;
  p.y = -determinant(l1.a, l1.c, l2.a, l2.c) / det;
  return p;
}

Line perpendicular_line(const Point &p, const Line &l) {
  Line pl = Line(l.b, -l.a, 0);
  pl.c = -(pl.a * p.x + pl.b * p.y);
  return pl.normed();
}

bool are_colinear(const Point &p1, const Point &p2, const Point &p3) {
  if (p1 == p2) return true;
  return line_by_two_points(p1, p2).contains(p3);
}

Point midpoint(const Point &p1, const Point &p2) { return Point((p1.x + p2.x) / 2., (p1.y + p2.y) / 2.); }

Circle circle_by_three_points();
Circle circle_by_center_and_point();

}  // nameespace tools
}  // namespace euclid_plane

using namespace euclid_plane;

class DumbStorage {
public:
  using Id = size_t;
  static constexpr Id kNoId = static_cast<Id>(-1);

  const Object* operator[](Id id) const { return pool_[id].get(); }

  template <typename T>
  Id add(T o) {
    Id id = find(o);
    if (id != kNoId) { return id; }
    add_new(o);
    return pool_.size() - 1;
  }

  template <typename T>
  T get(Id id) {
    Object* ptr = pool_[id].get();
    auto p = dynamic_cast<T*>(ptr);
    assert(p != nullptr);
    return *p;
  }

  // we need index here
  template <typename T>
  Id find(T o) const {
    auto type = o.type();
    for (Id i = 0; i < pool_.size(); ++i) {
      if (dumb_cmp(&o, pool_[i].get())) { return i; }
    }
    return kNoId;
  }

  std::vector<Id> get_ids_by_type(Type type) {
    std::vector<Id> ids;
    for (Id i = 0; i < pool_.size(); ++i) {
      if (pool_[i]->type() == type) { ids.push_back(i); }
    }
    return ids;
  }

  std::string to_json() {
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> w(buffer);

    w.StartObject();
    w.Key("commands");

    w.StartArray();
    for (auto& ptr : pool_) {
      std::string command;

      if (ptr->type() == Type::Line) {
        Line l = get_from_ptr<Line>(ptr.get());
        command = "(" + std::to_string(l.a) + ") * x + " +
            "(" + std::to_string(l.b) + ") * y + " +
            "(" + std::to_string(l.c) + ") = 0";
      }
      if (ptr->type() == Type::Point) {
        Point p = get_from_ptr<Point>(ptr.get());
        command = "Point({" + std::to_string(p.x) + ", " + std::to_string(p.y) + "})";
      }
      if (ptr->type() == Type::Circle) { assert(false); }

      w.String(command.c_str());
    }
    w.EndArray();

    w.EndObject();
    w.Flush();

    return std::string(buffer.GetString());
  }


private:
  template <typename T>
  T get_from_ptr(Object* o) {
    return *dynamic_cast<T*>(o);
  }

  template <typename T>
  void add_new(T o) {
    // do stuff
    pool_.push_back(std::make_unique<T>(o));
  }

  template<typename T>
  bool shit(Object* first, Object* second) const {
    auto f = dynamic_cast<T*>(first);
    auto s = dynamic_cast<T*>(second);
    if (f != nullptr && s != nullptr) {
      return *f == *s;
    }
    return false;
  }

  bool dumb_cmp(Object* first, Object* second) const {
    if (first->type() != second->type()) { return false; }
    if (shit<Point>(first, second)) { return true; }
    if (shit<Line>(first, second)) { return true; }
    if (shit<Circle>(first, second)) { return true; }
    return false;
  }

  std::vector<std::unique_ptr<Object>> pool_;

};


// enum class VertexTag {
//     Point,
//     Line,
//     Circle,
//     Last
// };

using VertexTag = Type;

enum class EdgeTag {
  Incidence = 0,
  Perpendicular,
  Parallel,
  Last
};



class Tool {
public:
  // virtual
  // arguments
  // some restrictions on arguments
  // new edges or vertices with edges

private:

};


class Graph {
public:
  using Id = size_t;
  struct Edge {
    Edge(Id v, Id u, EdgeTag tag) :
        v(v), u(u), tag(tag) {}
    Id v, u;
    EdgeTag tag;
    bool operator<(const Edge& other) const {
      return std::tie(v, u, tag) < std::tie(other.v, other.u, other.tag);
    }
  };
  using Edges = std::set<Edge>;

  std::pair<size_t, size_t> shape() const {
    return { tag_by_vertex_id_.size(), edges_.size() };
  }

  void add_vertex(Id v, VertexTag tag, int64_t comp = 0) {
    if (tag_by_vertex_id_.count(v)) { return; }
    complexity_[v] = comp;
    tag_by_vertex_id_[v] = tag;
  }

  void add_edge(Id v, Id u, EdgeTag tag) {
    if (v > u) { std::swap(v, u); }
    auto edge = Edge(v, u, tag);
    if (!edges_.count(edge)) {
      auto [it, _] = edges_.insert(edge);
      edges_by_vertex_id_[v].push_back(it);
      edges_by_vertex_id_[u].push_back(it);
    }
  }

  int64_t complexity(Id id) { return complexity_[id]; }

  std::vector<Id> get_ids(VertexTag tag, int64_t max_complexity) {
    std::vector<Id> ids;
    for (auto& [id, cmp] : complexity_) {
      if (cmp <= max_complexity && tag_by_vertex_id_[id] == tag) { ids.push_back(id); }
    }
    return ids;
  }


  std::string to_json() {
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> w(buffer);

    w.StartObject();
    w.Key("v");
    w.StartArray();
    for (auto& [id, tag] : tag_by_vertex_id_) {
      w.StartArray();
      w.Int64(id);
      w.Int64(static_cast<int64_t>(tag));
      w.Int64(complexity(id));
      w.EndArray();
    }
    w.EndArray();

    w.Key("e");
    w.StartArray();
    for (auto& e : edges_) {
      w.StartArray();
      w.Int64(e.v);
      w.Int64(e.u);
      w.Int64(static_cast<int64_t>(e.tag));
      w.EndArray();
    }
    w.EndArray();

    w.EndObject();
    w.Flush();

    return std::string(buffer.GetString());

    std::string ans;
    ans += "shape: " + std::to_string(tag_by_vertex_id_.size()) + ":" + std::to_string(edges_.size()) + "\n";
    ans += "V: ";
    for (auto& [k, v] : tag_by_vertex_id_) {
      ans += " { id: " + std::to_string(k) + " , tag: " + std::to_string(static_cast<size_t>(v)) + " comp: " + std::to_string(complexity_[k]) + " } ";
    }
    ans += "\nE: ";

    for (auto& e : edges_) {
      ans += " ( " + std::to_string(e.v) + ", " + std::to_string(e.u) + ", " + std::to_string(static_cast<size_t>(e.tag)) + " ) ";
    }
    ans += "\n";

    return ans;
  }

private:
  Edges edges_;
  std::map<Id, VertexTag> tag_by_vertex_id_;
  std::map<Id, std::vector<Edges::iterator>> edges_by_vertex_id_;
  std::map<Id, int64_t> complexity_;
};


enum class ActionTag {
  LineByTwoPoints = 0,
  PerpendicularLine,
  IntersectTwoLines,
  MidPoint,
  Last
};


struct Action {
  explicit Action(ActionTag tag) : tag(tag) {}
  Action() : Action(ActionTag::Last) {}
  virtual ~Action() {}
  virtual std::vector<Type> signature() = 0;
  virtual bool check(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const = 0;
  virtual std::vector<Graph::Id> operator()(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const = 0;
  ActionTag tag;
};

struct LineByTwoPoints : Action {
  LineByTwoPoints() : Action(ActionTag::LineByTwoPoints) {}
  virtual std::vector<Type> signature() { return { Type::Point, Type::Point }; }

  bool check(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const {
    if (arguments.size() != 2) { return false; }
    if (storage[arguments[0]]->type() != Type::Point) { return false; }
    if (storage[arguments[1]]->type() != Type::Point) { return false; }
    if (arguments[0] == arguments[1]) { return false; }

    return true;
  }

  std::vector<Graph::Id> operator()(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const {
    assert(check(graph, storage, arguments));

    const Object* o1 = storage[arguments[0]];
    auto ptr1 = dynamic_cast<const Point*>(o1);
    Point p1 = *ptr1;

    const Object* o2 = storage[arguments[1]];
    auto ptr2 = dynamic_cast<const Point*>(o2);
    Point p2 = *ptr2;

    Line l = tools::line_by_two_points(p1, p2);
    DumbStorage::Id id = storage.add(l);


    int64_t comp = 0;
    for (Id id : arguments) { comp = std::max(comp, graph.complexity(id)); }
    ++comp;

    graph.add_vertex(id, VertexTag::Line, comp);
    graph.add_edge(id, arguments[0], EdgeTag::Incidence);
    graph.add_edge(id, arguments[1], EdgeTag::Incidence);

    return { id };
  }
};

struct MidPoint : Action {
  MidPoint() : Action(ActionTag::MidPoint) {}
  virtual std::vector<Type> signature() { return { Type::Point, Type::Point }; }

  bool check(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const {
    if (arguments.size() != 2) { return false; }
    if (storage[arguments[0]]->type() != Type::Point) { return false; }
    if (storage[arguments[1]]->type() != Type::Point) { return false; }
    if (arguments[0] == arguments[1]) { return false; }

    return true;
  }

  std::vector<Graph::Id> operator()(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const {
    assert(check(graph, storage, arguments));

    const Object* o1 = storage[arguments[0]];
    auto ptr1 = dynamic_cast<const Point*>(o1);
    Point p1 = *ptr1;

    const Object* o2 = storage[arguments[1]];
    auto ptr2 = dynamic_cast<const Point*>(o2);
    Point p2 = *ptr2;

    Point p = tools::midpoint(p1, p2);
    DumbStorage::Id id = storage.add(p);

    int64_t comp = 0;
    for (Id id : arguments) { comp = std::max(comp, graph.complexity(id)); }
    ++comp;

    graph.add_vertex(id, VertexTag::Point, comp);
    // TODO(vk): EDGES?? EdgeTag HalfwayThere :D

    return { id };
  }
};

struct PerpendicularLine : Action {
  PerpendicularLine() : Action(ActionTag::PerpendicularLine) {}

  virtual std::vector<Type> signature() { return { Type::Point, Type::Line }; }

  bool check(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const {
    if (arguments.size() != 2) { return false; }
    if (storage[arguments[0]]->type() != Type::Point) { return false; }
    if (storage[arguments[1]]->type() != Type::Line) { return false; }
    return true;
  }

  std::vector<Graph::Id> operator()(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const {
    assert(check(graph, storage, arguments));

    const Object* o1 = storage[arguments[0]];
    auto ptr1 = dynamic_cast<const Point*>(o1);
    Point p1 = *ptr1;

    const Object* o2 = storage[arguments[1]];
    auto ptr2 = dynamic_cast<const Line*>(o2);
    Line l1 = *ptr2;

    Line l = tools::perpendicular_line(p1, l1);

    DumbStorage::Id id = storage.add(l);

    int64_t comp = 0;
    for (Id id : arguments) { comp = std::max(comp, graph.complexity(id)); }
    ++comp;

    graph.add_vertex(id, VertexTag::Line, comp);
    graph.add_edge(id, arguments[0], EdgeTag::Incidence);
    graph.add_edge(id, arguments[1], EdgeTag::Perpendicular);

    return { id };
  }
};


struct IntersectTwoLines: Action {
  IntersectTwoLines() : Action(ActionTag::IntersectTwoLines) {}

  virtual std::vector<Type> signature() { return { Type::Line, Type::Line }; }

  bool check(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const {
    if (arguments.size() != 2) { return false; }
    if (storage[arguments[0]]->type() != Type::Line) { return false; }
    if (storage[arguments[1]]->type() != Type::Line) { return false; }
    if (arguments[0] == arguments[1]) { return false; }
    return true;
  }

  std::vector<Graph::Id> operator()(Graph& graph, DumbStorage& storage, const std::vector<Graph::Id>& arguments) const {
    assert(check(graph, storage, arguments));

    const Object* o1 = storage[arguments[0]];
    auto ptr1 = dynamic_cast<const Line*>(o1);
    Line l1 = *ptr1;

    const Object* o2 = storage[arguments[1]];
    auto ptr2 = dynamic_cast<const Line*>(o2);
    Line l2 = *ptr2;

    if (tools::are_parallel(l1, l2)) {
      // TODO: add parallel relation
      graph.add_edge(arguments[0], arguments[1], EdgeTag::Parallel);
      return { };
    }

    Point p = tools::intersect_two_lines(l1, l2);


    DumbStorage::Id id = storage.add(p);

    int64_t comp = 0;
    for (Id id : arguments) { comp = std::max(comp, graph.complexity(id)); }
    ++comp;

    graph.add_vertex(id, VertexTag::Point, comp);
    graph.add_edge(id, arguments[0], EdgeTag::Incidence);
    graph.add_edge(id, arguments[1], EdgeTag::Incidence);

    return { id };
  }
};


class Backend {
public:
  using Id = DumbStorage::Id;
  Backend() {}
  void run() {
    init();
    // generate_noise();

    std::cout << graph_.to_json() << std::endl;
    std::cout << storage_.to_json() << std::endl;

    // std::ofstream f;
    // f.open("graph.json");
    // f << graph_.to_json();
    // f.close();

    // f.open("storage.json");
    // f << storage_.to_json();
    // f.close();

    while (true) {
      std::string raw_query;
      std::cin >> raw_query;

      rapidjson::Document json_query;
      json_query.Parse(raw_query.c_str());

      Graph::Id action_id = json_query["action_id"].GetInt();
      if (action_id == DumbStorage::kNoId) { break; }

      ActionTag action_tag = static_cast<ActionTag>(action_id);
      std::vector<Graph::Id> params;
      for (auto& v : json_query["params"].GetArray()) { params.push_back(v.GetInt()); }
      auto new_ids = action(action_tag, params);

      std::string json_state;
      json_state = graph_.to_json();

      // std::cout << json_state << std::endl;
      std::cout << new_ids << " " << graph_.to_json() << std::endl;
    }

  }

  void close() { std::cerr << "closed" << std::endl; }

  std::vector<Graph::Id> action(ActionTag tag, const std::vector<Graph::Id>& arguments) {
    // TODO: add log
    auto& a = *action_by_tag_[tag].get();
    if (!a.check(graph_, storage_, arguments)) { return {}; }
    auto ids = a(graph_, storage_, arguments);
    return ids;
  }

private:
  void generate_noise() {
    const int64_t max_complexity = 3;
    const int64_t iters = 0;
    for (size_t i = 0; i < iters; ++i) {
      std::map<Type, std::vector<Id>> ids_by_type;
      for (auto tag : { VertexTag::Point, VertexTag::Line }) {
        ids_by_type[tag] = graph_.get_ids(tag, max_complexity);
      }

      for (auto& [t, a] : action_by_tag_) {
        if (t == ActionTag::MidPoint) { continue; }
        auto signature = a->signature();
        std::vector<Id> args;
        bool ok = true;
        for (Type arg_t : signature) {
          auto& pool = ids_by_type[arg_t];
          if (pool.empty()) { ok = false; break; }
          int r = rnd() % pool.size();
          args.push_back(pool[r]);
        }
        if (ok) { action(t, args); }
      }

      // std::cerr << i << " " << graph_.shape() << std::endl;
    }

  }

  void init() {
    {
      add_action_to_pool<LineByTwoPoints>();
      add_action_to_pool<PerpendicularLine>();
      add_action_to_pool<IntersectTwoLines>();
      add_action_to_pool<MidPoint>();
    }

    std::vector<Point> ps = {{-1.9, 6.3}, {5, -3}, {-4, -3}};
    auto A = ps[0];
    auto B = ps[1];
    auto C = ps[2];

    storage_.add(A);
    storage_.add(B);
    storage_.add(C);

    graph_.add_vertex(0, VertexTag::Point, 0);
    graph_.add_vertex(1, VertexTag::Point, 0);
    graph_.add_vertex(2, VertexTag::Point, 0);

    bool add_h = true;
    if (add_h) {
      auto ab_id = action(ActionTag::LineByTwoPoints, {0, 1})[0];
      auto ac_id = action(ActionTag::LineByTwoPoints, {0, 2})[0];
      auto bc_id = action(ActionTag::LineByTwoPoints, {1, 2})[0];

      auto ha_id = action(ActionTag::PerpendicularLine, {0, bc_id})[0];
      auto hb_id = action(ActionTag::PerpendicularLine, {1, ac_id})[0];
      auto hc_id = action(ActionTag::PerpendicularLine, {2, ab_id})[0];

      auto Ha_id = action(ActionTag::IntersectTwoLines, {ha_id, bc_id})[0];
      auto Hb_id = action(ActionTag::IntersectTwoLines, {hb_id, ac_id})[0];
      auto Hc_id = action(ActionTag::IntersectTwoLines, {hc_id, ab_id})[0];

      auto h_id = action(ActionTag::IntersectTwoLines, {ha_id, hb_id})[0];

      auto ma_id = action(ActionTag::MidPoint, {1, 2})[0];
      auto l_id = action(ActionTag::LineByTwoPoints, {ma_id, h_id})[0];
      auto t_a_id = action(ActionTag::LineByTwoPoints, {Hb_id, Hc_id})[0];
      auto q_id = action(ActionTag::IntersectTwoLines, {t_a_id, bc_id})[0];

      auto aq_id = action(ActionTag::LineByTwoPoints, {q_id, 0})[0];
      auto l2_id = action(ActionTag::PerpendicularLine, {h_id, aq_id})[0];
      assert(l_id == l2_id);

      auto kc = action(ActionTag::PerpendicularLine, {2, ac_id})[0];

      auto d_id = action(ActionTag::IntersectTwoLines, {kc, l_id})[0];
      auto ad_id = action(ActionTag::LineByTwoPoints, {0, d_id})[0];

      // std::cout << graph_.to_json() << std::endl;
      // std::cout << graph_.shape() << std::endl;
    }

  }

  template<typename ActionType>
  void add_action_to_pool() {
    auto a = std::make_unique<ActionType>();
    auto tag = a->tag;
    action_by_tag_[tag] = std::move(a);
  }

  void run_tests() {
    // heights intersection
    {
      std::vector<Point> ps = {{0, -1},
                               {1, -1},
                               {2, 3}};
      auto A = ps[0];
      auto B = ps[1];
      auto C = ps[2];
      auto ab = tools::line_by_two_points(A, B);
      auto ac = tools::line_by_two_points(A, C);
      auto bc = tools::line_by_two_points(B, C);
      auto ha = tools::perpendicular_line(A, bc);
      auto hb = tools::perpendicular_line(B, ac);
      auto hc = tools::perpendicular_line(C, ab);
      auto H = tools::intersect_two_lines(ha, hb);
      assert(hc.contains(H));
    }
    // storage
    {

      std::vector<Point> ps = {{0, -1},
                               {1, -1},
                               {2, 3}};
      auto A = ps[0];
      auto B = ps[1];
      auto C = ps[2];
      auto ab = tools::line_by_two_points(A, B);
      auto ac = tools::line_by_two_points(A, C);
      auto bc = tools::line_by_two_points(B, C);
      auto ha = tools::perpendicular_line(A, bc);
      auto hb = tools::perpendicular_line(B, ac);
      auto hc = tools::perpendicular_line(C, ab);
      auto H = tools::intersect_two_lines(ha, hb);
      auto H1 = tools::intersect_two_lines(ha, hc);

      auto storage = DumbStorage();

      std::cout << storage.add(H) << std::endl;
      std::cout << storage.find(H1) << std::endl;
      std::cout << storage.add(A) << std::endl;
      std::cout << storage.add(bc) << std::endl;
      std::cout << storage.to_json() << std::endl;

    }
  }

  std::map<ActionTag, std::unique_ptr<Action>> action_by_tag_;

  DumbStorage storage_;
  Graph graph_;
};
