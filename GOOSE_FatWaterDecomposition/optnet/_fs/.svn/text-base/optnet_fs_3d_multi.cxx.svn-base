/*
 ==========================================================================
 |   
 |   $Id: optnet_fs_3d_multi.cxx 2137 2007-07-02 03:26:31Z kangli $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   Department of Electrical and Computer Engineering
 |   Carnegie Mellon University
 |   
 ==========================================================================
 |   This file is a part of the OptimalNet library.
 ==========================================================================
 | Copyright (c) 2003-2007 Kang Li <kangl@cmu.edu>. All rights reserved.
 | 
 | Version: MPL 1.1/GPL 2.0/LGPL 2.1
 | 
 | The contents of  this file are  subject to the  Mozilla Public License
 | Version  1.1 (the  "License"); you  may not  use this  file except  in
 | compliance with the License. You may obtain a copy of the License at
 | 
 | http://www.mozilla.org/MPL/
 | 
 | Software distributed under  the License is  distributed on an  "AS IS"
 | basis, WITHOUT WARRANTY  OF ANY KIND,  either express or  implied. See
 | the License for the specific language governing rights and limitations
 | under the License.
 | 
 | The Original Code is OptimalNet (optnet) Library code.
 | 
 | The  Initial  Developer of  the  Original Code  is  Kang Li.  Portions
 | created  by  the Initial  Developer  are Copyright  (C)  2003-2007 the
 | Initial Developer. All Rights Reserved.
 | 
 | Contributor(s): None
 | 
 | Alternatively, the contents of this  file may be used under  the terms
 | of either of the  GNU General Public License  Version 2 or later  (the
 | "GPL"), or the GNU Lesser General Public License Version 2.1 or  later
 | (the "LGPL"), in which case the provisions of the GPL or the LGPL  are
 | applicable instead of those  above. If you wish  to allow use of  your
 | version of this  file only under  the terms of  either the GPL  or the
 | LGPL, and not to allow others  to use your version of this  file under
 | the  terms  of  the  MPL,  indicate  your  decision  by  deleting  the
 | provisions above and replace them with the notice and other provisions
 | required by the GPL or the  LGPL. If you do not delete  the provisions
 | above, a recipient may use your  version of this file under the  terms
 | of any one of the MPL, the GPL or the LGPL.
 ==========================================================================
 */

#ifndef ___OPTNET_FS_3D_MULTI_CXX___
#   define ___OPTNET_FS_3D_MULTI_CXX___

#   include <optnet/_base/except.hxx>
#   include <optnet/_fs/optnet_fs_3d_multi.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <deque>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::optnet_fs_3d_multi() :
    m_pcost(0)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::create(
                                        const cost_array_base_type& cost
                                        )
{
    if (!m_graph.create(cost.size_0(), 
                        cost.size_1(), 
                        cost.size_2(),
                        cost.size_3()
                        )) {
        throw_exception(std::runtime_error(
            "optnet_fs_3d_multi::create: Could not create graph."
        ));
    }

    m_intra.resize(m_graph.size_3());

    // Set up default parameters.
    for (size_type i = 0; i < m_graph.size_3(); ++i) {
    m_intra[i].tonode[0].clear();
    m_intra[i].tonode[1].clear();
        m_intra[i].dropped   = false;
        m_intra[i].climbed   = false;
        m_intra[i].circle[0] = false;
        m_intra[i].circle[1] = false;
        m_intra[i].margin[0] = 0;
        m_intra[i].margin[1] = 0;
        m_intra[i].smooth[0] = 1;
        m_intra[i].smooth[1] = 1;
    }

    m_inter.clear();

    // Temporarily save a pointer to the cost vector.
    m_pcost = &cost;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::set_params(size_type k,
                                                 int       smooth0,
                                                 int       smooth1,
                                                 bool      circle0,
                                                 bool      circle1
                                                 )
{
    // Check arguments.
    if (k >= m_intra.size() ||  smooth0 < 0 || smooth1 < 0) {
        throw_exception(
            std::invalid_argument(
            "optnet_fs_3d_multi::set_params: Invalid argument."
        ));
    }

    m_intra[k].smooth[0] = smooth0;
    m_intra[k].smooth[1] = smooth1;
    m_intra[k].circle[0] = circle0;
    m_intra[k].circle[1] = circle1;

}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::set_relation(size_type k0,
                                                   size_type k1,
                                                   int       r01,
                                                   int       r10
                                                   )
{
    // Check arguments.
    if ((k0 >= m_graph.size_3()) || (k1 >= m_graph.size_3()) ||
        (r01 >= 0 && r10 >= 0)) {
        throw_exception(
            std::invalid_argument(
            "optnet_fs_3d_multi::set_relation: Invalid argument."
        ));
    }
    
    // The relations poses constraints that some nodes definitely cannot
    //  be on the final surfaces. 
    //
    // Note that:
    //  if r01 < 0 || r10 < 0 but not both : non-crossing case
    //  if r01 < 0 && r10 < 0              : crossing case
    //  if r01 > 0 && r10 > 0              : invalid
    //
    if (r01 * r10 <= 0) {
        
        // Negative number specifies maximum distance, possitive
        // number specifies minimum distance.
        if (r01 + r10 > 0) {
            throw_exception(
                std::invalid_argument(
                "optnet_fs_3d_multi::set_relation: Invalid surface interrelation."
            ));
        }

        // For lower margin.
        m_intra[k0].tonode[0].push_back(std::pair<size_t, int>(k1, r01));
        m_intra[k1].tonode[0].push_back(std::pair<size_t, int>(k0, r10));

        // For upper margin.
        m_intra[k1].tonode[1].push_back(std::pair<size_t, int>(k0, r01));
        m_intra[k0].tonode[1].push_back(std::pair<size_t, int>(k1, r10));

        if (r01 > 0) {
            m_intra[k1].climbed = true;
            m_intra[k0].dropped = true;
        }
        else if (r10 > 0) {
            m_intra[k0].climbed = true;
            m_intra[k1].dropped = true;
        }
    }

    // Create a new relation and 
    //  append it to the relation vector.
    _Inter relation;

    relation.r[0] = r01;
    relation.r[1] = r10;
    relation.k[0] = k0;
    relation.k[1] = k1;
    
    m_inter.push_back(relation);

}

#   ifdef __OPTNET_SUPPORT_ROI__

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::set_roi(const roi_base_type& roi)
{
    if (!m_graph.set_roi(roi)) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_fs_3d_multi::set_roi: The ROI size must match the graph size."
        ));
    }
}

#   endif // __OPTNET_SUPPORT_ROI__

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::solve(net_base_type& net, 
                                            capacity_type* pflow
                                            )
{
    size_type       i0, i1, i2, i3;
    capacity_type   flow;

    if (net.size_0() != m_graph.size_0() || 
        net.size_1() != m_graph.size_1() ||
        net.size_2() != m_graph.size_3()
        ) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_fs_3d_multi::solve: The output surface size must match the graph size."
        ));
    }

    // Compute the lower and upper margin of the graph nodes in
    // each 3-D subgraph.
    get_bounds_of_subgraphs();

    // Assign the cost of graph nodes based on the input cost
    // vector. We also perform the "translation operation"
    // here to guaranttee a non-empty solution.
    transform_costs();

    // Build the arcs of the graphs.
    build_arcs();

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();

    // Recover optimal surface(s) from source set.
    for (i3 = 0; i3 < m_graph.size_3(); ++i3) {
        for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
            for (i0 = 0; i0 < m_graph.size_0(); ++i0) {

                size_type lowest;
                
#   ifdef __OPTNET_SUPPORT_ROI__
                if (m_graph.has_roi()) {
                    lowest = m_graph.roi_node(i0, i1, i3).lower;
                    if (m_intra[i3].margin[0] > (int)lowest)
                        lowest = (size_type)m_intra[i3].margin[0];
                } else
#   endif
                    lowest = m_intra[i3].margin[0];

                for (i2 = lowest;
                     i2 < m_graph.size_2() - m_intra[i3].margin[1];
                     ++i2) {
                    
                    // Find upper envelope.
                    if (!m_graph.in_source_set(i0, i1, i2, i3))
                        break;
                }
                if (i2 != lowest)
                    net(i0, i1, i3) = (net_base_type::value_type)(i2 - 1);
                else
                    net(i0, i1, i3) = (net_base_type::value_type)(-1);
            } // for i0
        } // for i1
    } // for i3

    if (0 != pflow)
        *pflow = flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::get_bounds_of_subgraphs()
{
    bool                stop[2];
    int                 r, new_bound;
    size_type           i, k, n = m_intra.size();
    std::deque<_Intra*> Q[2];

    //
    // Hopefully this is also a correct way of detecting illegal surface
    // relation specifications. 
    //

    //
    // STEP 1: Initialize.
    //
    for (i = 0; i < n; ++i) {

        _Intra* p = &m_intra[i];
        
        if (p->climbed) p->visits[0] = 0;
        else { // else 1
            Q[0].push_back(p);
            p->visits[0] = 1;
        } // else 1

        if (p->dropped) p->visits[1] = 0;
        else { // else 2
            Q[1].push_back(p);
            p->visits[1] = 1;
        } // else 2

    } // for

    //
    // STEP 2: Compute the margin of the graphs.
    //
    if (Q[0].empty() || Q[1].empty()) {
        // The specified relations are invalid.
        throw_exception(
            std::invalid_argument(
            "optnet_fs_3d_multi::get_bounds_of_subgraphs: Invalid surface interrelation."
        ));
    }
    else {
        stop[0] = false;
        stop[1] = false;
        
        //
        // Lower margin.
        do {
            _Intra* p = Q[0].front();
            Q[0].pop_front();

            for (i = 0; i < p->tonode[0].size(); ++i) {
                k = p->tonode[0][i].first;
                r = p->tonode[0][i].second;
                new_bound = p->margin[0] + r;
                if (m_intra[k].margin[0] < new_bound) {
                    if (++m_intra[k].visits[0] < (int)n) {
                        m_intra[k].margin[0] = new_bound;
                        Q[0].push_back(&m_intra[k]);
                    }
                    else {
                        stop[0] = true;
                        break;
                    } // else
                } // if
            } // for

        } // do
        while (!Q[0].empty() && !stop[0]);

        //
        // Upper margin.
        do {
            _Intra* p = Q[1].front();
            Q[1].pop_front();

            for (i = 0; i < p->tonode[1].size(); ++i) {
                k = p->tonode[1][i].first;
                r = p->tonode[1][i].second;
                new_bound = p->margin[1] + r;
                if (m_intra[k].margin[1] < new_bound) {
                    if (++m_intra[k].visits[1] < (int)n) {
                        m_intra[k].margin[1] = new_bound;
                        Q[1].push_back(&m_intra[k]);
                    }
                    else {
                        stop[1] = true;
                        break;
                    } // else
                } // if
            } // for

        } // do
        while (!Q[1].empty() && !stop[1]);

        if (stop[0] || stop[1]) {
            // There is a loop in the relations specified.
            throw_exception(
                std::invalid_argument(
                "optnet_fs_3d_multi::get_bounds_of_subgraphs: Invalid surface interrelation."
            ));
        } // if
    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::transform_costs()
{
    assert(m_pcost != 0);

    size_type            i0, i1, i2, i3;

    const size_type& s0 = m_pcost->size_0();
    const size_type& s1 = m_pcost->size_1();
    const size_type& s2 = m_pcost->size_2();
    const size_type& s3 = m_pcost->size_3();

    m_graph.set_initial_flow(0);


#   ifdef __OPTNET_SUPPORT_ROI__

    if (m_graph.has_roi()) {

        for (i3 = 0; i3 < s3; ++i3) {
            // Nodes below margin[0] and above s2-margin[1] will
            // be ignored.
            for (i2 = m_intra[i3].margin[0];
                i2 < s2 - m_intra[i3].margin[1];
                ++i2) {

                for (i1 = 0; i1 < s1; ++i1) {
                    for (i0 = 0; i0 < s0; ++i0) {

                        if (m_graph.in_roi(i0, i1, i2, i3)) {
                            
                            if (((int)i2 == m_intra[i3].margin[0]) ||
                                (m_graph.is_lowest(i0, i1, i2, i3))
                                ) {
                                // The node is the lowest in the column, make its cost -1.
                                m_graph.add_st_arc(1, 0, i0, i1, i2, i3);
                            }
                            else {

                                capacity_type cap
                                    = (capacity_type)(*m_pcost)(i0, i1, i2, i3)
                                    - (capacity_type)(*m_pcost)(i0, i1, i2 - 1, i3);

                                if (cap >= 0) // non-negative -> connect to t
                                    m_graph.add_st_arc(0, +cap, i0, i1, i2, i3);
                                else          // negative     -> connect to s
                                    m_graph.add_st_arc(-cap, 0, i0, i1, i2, i3);
                        
                            }

                        } // if (m_graph.in_roi(...

                    } // for i0
                } // for i1

            } // for i2
        } // for i3     
        
    }
    else { // else 1

#   endif // __OPTNET_SUPPORT_ROI__

        // Construct the s-t graph "G_st".
        for (i3 = 0; i3 < s3; ++i3) {
            // Nodes below margin[0] and above s2-margin[1] will
            // be ignored.
            for (i2 = m_intra[i3].margin[0] + 1;
                i2 < s2 - m_intra[i3].margin[1];
                ++i2) {
                for (i1 = 0; i1 < s1; ++i1) {
                    for (i0 = 0; i0 < s0; ++i0) {

                        capacity_type cap
                            = (capacity_type)(*m_pcost)(i0, i1, i2, i3)
                            - (capacity_type)(*m_pcost)(i0, i1, i2 - 1, i3);

                        if (cap >= 0) // non-negative -> connect to t
                            m_graph.add_st_arc(0, +cap, i0, i1, i2, i3);
                        else          // negative     -> connect to s
                            m_graph.add_st_arc(-cap, 0, i0, i1, i2, i3);

                    } // for i0
                } // for i1
            } // for i2
        } // for i3

        for (i3 = 0; i3 < s3; ++i3) {
            i2 = m_intra[i3].margin[0]; // lowest margin
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {
                    // cost = -1
                    m_graph.add_st_arc(1, 0, i0, i1, i2, i3);
                } // for i0
            } // for i1
        } // for i3


#   ifdef __OPTNET_SUPPORT_ROI__
    } // else 1
#   endif // __OPTNET_SUPPORT_ROI__

}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d_multi<_Cost, _Cap, _Tg>::build_arcs()
{
    int ii, i0, i1, i2, i3, i2a, i2b;
    int s0, s1, s2, s3;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();
    s3 = (int)m_graph.size_3();

    // Clear all constructed arcs.
    m_graph.clear_arcs();

    //
    // Construct intra-surface arcs here.
    //

    for (i3 = 0; i3 < s3; ++i3) {

        const bool& circle0 = m_intra[i3].circle[0];
        const bool& circle1 = m_intra[i3].circle[1];
        const int & bounds0 = m_intra[i3].margin[0];
        const int & bounds1 = m_intra[i3].margin[1];
        const int & smooth0 = m_intra[i3].smooth[0];
        const int & smooth1 = m_intra[i3].smooth[1];

        // -- Intra-column (vertical) arcs.
        for (i2 = s2 - bounds1 - 1;
             i2 > bounds0;
             --i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1,          i2 - 1,          i3);
                } // for i0 
            } // for i1
        } // for i2

        // MKH -- only add arcs if more than one layer in this direction
        if (s0 > 1)
        {
        // -- Inter-column arcs (dir-0).
        for (i2 = s2 - bounds1 - 1;
             i2 > smooth0 + bounds0;
             --i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 1; i0 < s0 - 1; ++i0) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 - 1,      i1,          i2 - smooth0,    i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 + 1,      i1,          i2 - smooth0,    i3);
                } {
                    m_graph.add_arc(0,           i1,          i2,          i3,          1,           i1,          i2 - smooth0,    i3);
                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          s0 - 2,      i1,          i2 - smooth0,    i3);
                } // for i0
            } // for i1
        } // for i2
        }
        /////

        // MKH -- only add arcs if more than one layer in this direction
        if (s1 > 1)
        {
        // -- Inter-column arcs (dir-1).
        for (i2 = s2 - bounds1 - 1;
             i2 > smooth1 + bounds0;
             --i2) {
            for (i0 = 0; i0 < s0; ++i0) {
                for (i1 = 1; i1 < s1 - 1; ++i1) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 - 1,      i2 - smooth1,    i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 + 1,      i2 - smooth1,    i3);
                } {
                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          1,           i2 - smooth1,    i3);
                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          s1 - 2,      i2 - smooth1,    i3);
                } // for i0
            } // for i1
        } // for i2
        }
        /////

        // -- The base-set
        { 
            // MKH -- only add arcs if more than one layer in this direction
            if (s0 > 1)
            {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 1; i0 < s0 - 1; ++i0) {

#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(i0    , i1, i3).lower;
                        i2a = (int)m_graph.roi_node(i0 - 1, i1, i3).lower;
                        i2b = (int)m_graph.roi_node(i0 + 1, i1, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2a) i2a = bounds0;
                        if (bounds0 > i2b) i2b = bounds0;
                    } else
#   endif
                    i2  = i2a = i2b = bounds0;

                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 - 1,      i1,          i2a,             i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 + 1,      i1,          i2b,             i3);
                } {
                    // i0 = 0
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(     0, i1, i3).lower;
                        i2b = (int)m_graph.roi_node(     1, i1, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2b) i2b = bounds0;
                    } else
#   endif
                    i2  = i2b = bounds0;

                    m_graph.add_arc(0,           i1,          i2,          i3,          1,           i1,          i2b,             i3);

                    // i0 = s0 - 1
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(s0 - 1, i1, i3).lower;
                        i2a = (int)m_graph.roi_node(s0 - 2, i1, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2a) i2a = bounds0;
                    } else
#   endif
                        i2  = i2a = bounds0;

                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          s0 - 2,      i1,          i2a,             i3);
                } // for i0
            } // for i1
            }
            /////
        }

        {
            // MKH -- only add arcs if more than one layer in this direction
            if (s1 > 1)
            {
            for (i0 = 0; i0 < s0; ++i0) {
                for (i1 = 1; i1 < s1 - 1; ++i1) {
                    
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(i0, i1    , i3).lower;
                        i2a = (int)m_graph.roi_node(i0, i1 - 1, i3).lower;
                        i2b = (int)m_graph.roi_node(i0, i1 + 1, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2a) i2a = bounds0;
                        if (bounds0 > i2b) i2b = bounds0;
                    } else
#   endif
                        i2  = i2a = i2b = bounds0;

                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 - 1,      i2a,             i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 + 1,      i2b,             i3);
                } {
                    // i1 = 0
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(i0,      0, i3).lower;
                        i2b = (int)m_graph.roi_node(i0,      1, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2b) i2b = bounds0;
                    } else
#   endif
                        i2  = i2b = bounds0;

                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          1,           i2b,             i3);

                    // i1 = s1 - 1
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(i0, s1 - 1, i3).lower;
                        i2a = (int)m_graph.roi_node(i0, s1 - 2, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2a) i2a = bounds0;
                    } else
#   endif
                        i2  = i2a = bounds0;

                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          s1 - 2,      i2a,             i3);
                } // for i0
            } // for i1
            }
            /////
        }

        // -- Circular graph connections.
        if (circle0) {
            for (i2 = s2 - bounds1 - 1;
                 i2 > smooth0 + bounds0;
                 --i2) {
                for (i1 = 0; i1 < s1; ++i1) {
                    m_graph.add_arc(0,           i1,          i2,          i3,          s0 - 1,      i1,          i2 - smooth0,    i3);
                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          0,           i1,          i2 - smooth0,    i3);
                }
            } {
                for (i1 = 0; i1 < s1; ++i1) {
                    // i0 = 0
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(0     , i1, i3).lower;
                        i2a = (int)m_graph.roi_node(s0 - 1, i1, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2a) i2a = bounds0;
                    } else
#   endif
                        i2  = i2a = bounds0;

                    m_graph.add_arc(0,           i1,          i2,          i3,          s0 - 1,      i1,          i2a,             i3);
 
                    // i0 = s0 - 1
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(s0 - 1, i1, i3).lower;
                        i2b = (int)m_graph.roi_node(0     , i1, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2b) i2b = bounds0;
                    } else
#   endif
                        i2  = i2b = bounds0;

                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          0,           i1,          i2b,             i3);
                }
            }
        } // if

        if (circle1) {
            for (i2 = s2 - bounds1 - 1;
                 i2 > smooth1 + bounds0;
                 --i2) {
                for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          s1 - 1,      i2 - smooth1,    i3);
                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          0,           i2 - smooth1,    i3);
                }                                                          
            } {                                      
                for (i0 = 0; i0 < s0; ++i0) {                              
                    // i1 = 0
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(i0, 0     , i3).lower;
                        i2a = (int)m_graph.roi_node(i0, s1 - 1, i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2a) i2a = bounds0;
                    } else
#   endif
                        i2  = i2a = bounds0;

                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          s1 - 1,      i2a,             i3);

                    // i1 = s1 - 1
#   ifdef __OPTNET_SUPPORT_ROI__
                    if (m_graph.has_roi()) {
                        i2  = (int)m_graph.roi_node(i0, s1 - 1, i3).lower;
                        i2b = (int)m_graph.roi_node(i0, 0     , i3).lower;

                        if (bounds0 > i2 ) i2  = bounds0;
                        if (bounds0 > i2b) i2b = bounds0;
                    } else
#   endif
                        i2  = i2b = bounds0;

                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          0,           i2b,            i3);
                }
            }
        } // if

    } // for i3
    
    // Construct inter-surface arcs here.
    for (i3 = 0; i3 < (int)m_inter.size(); ++i3) {

        const size_type& k0 = m_inter[i3].k[0];
        const size_type& k1 = m_inter[i3].k[1];
        const int&       r0 = m_inter[i3].r[0];
        const int&       r1 = m_inter[i3].r[1];

        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {

                // k0 --> k1
                for (i2 = m_intra[k0].margin[0];
                     i2 < s2 - m_intra[k0].margin[1];
                     ++i2) {
                    
                    ii = i2 + r0;
                    if (ii <= m_intra[k1].margin[0] || 
                        ii >= s2 - m_intra[k1].margin[1]) continue;
                    m_graph.add_arc(i0, i1, i2, k0, i0, i1, ii, k1);
                }

                // k1 --> k0
                for (i2 = m_intra[k1].margin[0];
                     i2 < s2 - m_intra[k1].margin[1];
                     ++i2) {
                    
                    ii = i2 + r1;
                    if (ii <= m_intra[k0].margin[0] || 
                        ii >= s2 - m_intra[k0].margin[1]) continue;
                    m_graph.add_arc(i0, i1, i2, k1, i0, i1, ii, k0);
                }

            } // for i0
        } // for i1
    } // for i3

    // Interconnect "zero-plane"s.
    for (i3 = 1; i3 < s3; ++i3) {

        m_graph.add_arc(0, 0, m_intra[i3    ].margin[0], i3,     
                        0, 0, m_intra[i3 - 1].margin[0], i3 - 1
                        );
        m_graph.add_arc(0, 0, m_intra[i3 - 1].margin[0], i3 - 1,
                        0, 0, m_intra[i3    ].margin[0], i3
                        );
    }

}


} // namespace

#endif
