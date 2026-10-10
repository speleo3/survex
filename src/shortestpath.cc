//
//  Path finding.
//
//  Copyright (C) 2021 Thomas Holder
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "shortestpath.h"

#include "model.h"
#include "vector3.h"

#include <deque>
#include <map>
#include <vector>
#include <algorithm>

namespace svx {

using ConnectedStationKey = Vector3;

template <typename T>
std::vector<T const*> const& const_item_cast(std::vector<T*> const& v) {
    return reinterpret_cast<std::vector<T const*> const&>(v);
}

/**
 * Station with an adjacency list of connected stations.
 */
class ConnectedStation {
    using Connected_t = std::vector<ConnectedStation *>;
    using Connected_const_t = std::vector<ConnectedStation const*> const;

    ConnectedStationKey const& m_point;
    Connected_t m_connected;
    unsigned m_disjunct_set_index = 0;

  public:
    ConnectedStation(ConnectedStationKey const& p) : m_point(p) {}

    Connected_const_t & connected() const {
        return const_item_cast(m_connected);
    }

    Vector3 const& point() const { return m_point; }

    void connect(ConnectedStation& other) {
        assert(m_disjunct_set_index == 0);
        m_connected.push_back(&other);
        other.m_connected.push_back(this);
    }

    double distance(ConnectedStation const& other) const {
        return (m_point - other.m_point).magnitude();
    }

    unsigned disjunct_set_index() const { return m_disjunct_set_index; }

    bool add_to_disjunct_set(unsigned set_num) {
        assert(set_num > 0);

        if (m_disjunct_set_index != 0) {
            return false;
        }

        m_disjunct_set_index = set_num;

        for (auto* y : m_connected) {
            y->add_to_disjunct_set(set_num);
        }

        return true;
    }
};

/**
 * Set of stations. Stations in the set are unique by XYZ position.
 */
class ConnectedStationSet {
    std::map<ConnectedStationKey, ConnectedStation> m_data;

  public:
    ConnectedStation& emplace(ConnectedStationKey const& point) {
        return m_data.try_emplace(point, point).first->second;
    }

    ConnectedStation const* find(ConnectedStationKey const& point) const {
        auto const it = m_data.find(point);
        return it == m_data.end() ? nullptr : &it->second;
    }

    ConnectedStationSet(Model const& model) {
        for (int f = 0; f != 8; ++f) {
            for (auto it = model.traverses_begin(f, nullptr),
                      it_end = model.traverses_end(f);
                 it != it_end; ++it) {
                ConnectedStation* prev_station = nullptr;

                for (auto& point : *it) {
                    auto& station = this->emplace(point);

                    if (prev_station) {
                        prev_station->connect(station);
                    }

                    prev_station = &station;
                }
            }
        }

        unsigned current_set = 1;
        for (auto& item : m_data) {
            if (item.second.add_to_disjunct_set(current_set)) {
                ++current_set;
            }
        }
    }
};

/**
 * Shortest Path between two stations.
 *
 * Returns tuple of the length of the path (or -1 if no path is found) and
 * list of stations along the path.
 *
 * Uses A* algorithm, adapted from Wikipedia:
 * http://en.wikipedia.org/w/index.php?title=A*_search_algorithm&oldid=289896415
 */
static std::pair<double, std::vector<ConnectedStation const*>>
shortestpath(ConnectedStation const& self, ConnectedStation const& other) {
    if (self.disjunct_set_index() != other.disjunct_set_index()) {
        return {-1, {}};
    }

    auto const* const self_ptr = &self;
    auto const* const other_ptr = &other;

    auto came_from = std::map<ConnectedStation const*,
                              ConnectedStation const*>(); // Map for backtrace

    auto closedset = std::set<ConnectedStation const*>{}; // The set of nodes
                                                          // already evaluated.
    auto openset =
        std::deque<ConnectedStation const*>{self_ptr}; // List sorted by f_score
    auto g_score = std::map<ConnectedStation const*, double>{
        {self_ptr, 0}}; // Distance from self along optimal path.
    auto h_score = std::map<ConnectedStation const*, double>{
        {self_ptr,
         self.distance(other)}}; // Estimated lower bound from y to other
    auto f_score = std::map<ConnectedStation const*, double>{
        {self_ptr, h_score[self_ptr]}}; // Estimated total distance from self
                                        // to other through y.
    while (!openset.empty()) {
        auto const* const x = openset.front();
        openset.pop_front();

        if (x == other_ptr) {
            auto path = std::vector<ConnectedStation const*>{x};

            for (auto y = x; (y = came_from[y]);) {
                path.push_back(y);
            }

            std::reverse(path.begin(), path.end());

            return {g_score[other_ptr], path};
        }

        closedset.insert(x);

        for (auto const* const y : x->connected()) {
            if (closedset.count(y)) {
                continue;
            }

            auto tentative_g_score = g_score[x] + x->distance(*y);
            bool tentative_is_better = false;

            auto it = std::find(openset.begin(), openset.end(), y);
            if (it == openset.end()) {
                h_score[y] = y->distance(other);
                tentative_is_better = true;
            } else if (tentative_g_score < g_score[y]) {
                openset.erase(it);
                tentative_is_better = true;
            }

            if (tentative_is_better) {
                came_from[y] = x;
                g_score[y] = tentative_g_score;
                f_score[y] = tentative_g_score + h_score[y];

                for (it = openset.begin(); it != openset.end(); ++it) {
                    if (f_score[y] < f_score[*it]) {
                        break;
                    }
                }

                openset.insert(it, y);
            }
        }
    }

    return {-1, {}};
}

/**
 * Shortest Path between two stations.
 */
std::pair<double, std::vector<Vector3>>
shortestpath(Model const& model, Vector3 const& start, Vector3 const& end) {
    auto const stations = svx::ConnectedStationSet(model);
    auto* station_start = stations.find(start);
    auto* station_end = stations.find(end);

    if (!station_start || !station_end) {
        return {-1, {}};
    }

    auto [dist, path] = svx::shortestpath(*station_start, *station_end);

    std::vector<Vector3> vecpath;

    for (auto const& station : path) {
        vecpath.push_back(station->point());
    }

    return {dist, vecpath};
}

} // namespace svx
