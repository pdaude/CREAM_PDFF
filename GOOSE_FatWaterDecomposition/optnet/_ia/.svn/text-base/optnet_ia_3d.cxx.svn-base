/*
 ==========================================================================
 |   
 |   $Id: optnet_ia_3d.cxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___OPTNET_IA_3D_CXX___
#   define ___OPTNET_IA_3D_CXX___

#   include <optnet/_base/except.hxx>
#   include <optnet/_ia/optnet_ia_3d.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_ia_3d<_Cost, _Cap, _Tg>::optnet_ia_3d() :
    m_pcost(0)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_ia_3d<_Cost, _Cap, _Tg>::create(
                                    const cost_array_base_type& cost
                                    )
{
    if (!m_graph.create(cost.size_0(), 
                        cost.size_1(), 
                        cost.size_2()
        )) {
        throw_exception(std::runtime_error(
            "optnet_ia_3d::create: Could not create graph."
            ));
    }

    // Temporarily save a pointer to the cost vector.
    m_pcost = &cost;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_ia_3d<_Cost, _Cap, _Tg>::set_params(int  smooth0,
                                           int  smooth1,
                                           bool circle0,
                                           bool circle1
                                           )
{
    if (!m_graph.set_params(smooth0, smooth1, circle0, circle1)) {
        throw_exception(
            std::invalid_argument(
            "optnet_ia_3d::set_params: Invalid argument."
        ));
    } // if
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_ia_3d<_Cost, _Cap, _Tg>::set_roi(const roi_base_type& roi)
{
    if (!m_graph.set_roi(roi)) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_ia_3d::set_roi: The ROI size must match the graph size."
        ));
    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_ia_3d<_Cost, _Cap, _Tg>::solve(net_base_type& net,
                                      capacity_type* pflow)
{
    size_type       i0, i1, i2;
    capacity_type   flow;

    if (net.size_0() != m_graph.size_0() || 
        net.size_1() != m_graph.size_1()
        ) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_ia_3d::solve: The output surface size must match the graph size."
            ));
    }

    // Assign the cost of graph nodes based on the input cost
    // vector. We also perform the "translation operation"
    // here to guaranttee a non-empty solution.
    transform_costs();

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();

    // Recover optimal surface(s) from source set.
    for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
        for (i0 = 0; i0 < m_graph.size_0(); ++i0) {

            size_type lowest = m_graph.lowest(i0, i1);

            for (i2 = lowest; i2 < m_graph.size_2(); ++i2) {
                // Find upper envelope.
                if (!m_graph.in_source_set(i0, i1, i2))
                    break;
            }

            if (i2 != lowest)
                net(i0, i1, 0) = (int)(i2 - 1);
            else
                net(i0, i1, 0) = -1;
        }
    }

    if (0 != pflow)
        *pflow = flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_ia_3d<_Cost, _Cap, _Tg>::transform_costs()
{
    assert(m_pcost != 0);

    size_type            i0, i1, i2;

    const size_type& s0 = m_pcost->size_0();
    const size_type& s1 = m_pcost->size_1();
    const size_type& s2 = m_pcost->size_2();

    m_graph.set_initial_flow(0);

    //
    // Construct the s-t graph "G_st".
    //
    if (m_graph.has_roi()) {


        for (i2 = 0; i2 < s2; ++i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {

                    if (m_graph.in_roi(i0, i1, i2)) {

                        if (m_graph.is_lowest(i0, i1, i2)) {
                            // The node is the lowest in the column, make its cost -1.
                            m_graph.add_st_arc(1, 0, i0, i1, i2);
                        }
                        else {
                            capacity_type cap = (capacity_type)(*m_pcost)(i0, i1, i2)
                                              - (capacity_type)(*m_pcost)(i0, i1, i2 - 1);

                            if (cap >= 0) // non-negative -> connect to t
                                m_graph.add_st_arc(0, +cap, i0, i1, i2);
                            else          // negative     -> connect to s
                                m_graph.add_st_arc(-cap, 0, i0, i1, i2);
                        }

                    } // if (m_graph.in_roi(i0, i1, i2))

                } // for i0
            } // for i1
        } // for i2

    }
    else {

        for (i2 = 1; i2 < s2; ++i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {

                    capacity_type cap = (capacity_type)(*m_pcost)(i0, i1, i2)
                                      - (capacity_type)(*m_pcost)(i0, i1, i2 - 1);

                    if (cap >= 0) // non-negative -> connect to t
                        m_graph.add_st_arc(0, +cap, i0, i1, i2);
                    else          // negative     -> connect to s
                        m_graph.add_st_arc(-cap, 0, i0, i1, i2);
                    
                } // for i0
            } // for i1
        } // for i2

        // s2 = 0
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                m_graph.add_st_arc(1, 0, i0, i1, 0);
            } // for i0
        } // for i1

    } // else
}


} // namespace

#endif
