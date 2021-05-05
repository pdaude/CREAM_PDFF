/*
 ==========================================================================
 |   
 |   $Id: optnet_fs_3d.cxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___OPTNET_FS_3D_CXX___
#   define ___OPTNET_FS_3D_CXX___

#   include <optnet/_base/except.hxx>
#   include <optnet/_fs/optnet_fs_3d.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_fs_3d<_Cost, _Cap, _Tg>::optnet_fs_3d() :
    m_pcost(0)
{
    // Initialize smoothness and circleity parameters.
    m_smooth[0] = 1;
    m_smooth[1] = 1;
    m_circle[0] = false;
    m_circle[1] = false;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d<_Cost, _Cap, _Tg>::create(
                                    const cost_array_base_type& cost
                                    )
{
    if (!m_graph.create(cost.size_0(), 
                        cost.size_1(), 
                        cost.size_2()
                        )) {
        throw_exception(std::runtime_error(
            "optnet_fs_3d::create: Could not create graph."
            ));
    }

    // Save a pointer to the cost vector.
    m_pcost = &cost;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d<_Cost, _Cap, _Tg>::set_params(int  smooth0,
                                           int  smooth1,
                                           bool circle0,
                                           bool circle1
                                           )
{
    if (smooth0 < 0 || smooth1 < 0) {
        throw_exception(std::invalid_argument(
            "optnet_fs_3d::set_params: Invalid argument."
            ));
    }

    m_smooth[0] = smooth0;
    m_smooth[1] = smooth1;
    m_circle[0] = circle0;
    m_circle[1] = circle1;
}

#   ifdef __OPTNET_SUPPORT_ROI__

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d<_Cost, _Cap, _Tg>::set_roi(const roi_base_type& roi)
{
    if (!m_graph.set_roi(roi)) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_fs_3d::set_roi: The ROI size must match the graph size."
        ));
    }
}

#   endif // __OPTNET_SUPPORT_ROI__

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d<_Cost, _Cap, _Tg>::solve(net_base_type& net,
                                      capacity_type* pflow
                                      )
{
    size_type       i0, i1, i2;
    capacity_type   flow;

    if (net.size_0() != m_graph.size_0() || 
        net.size_1() != m_graph.size_1()
        ) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_fs_3d::solve: The net size must match the graph size."
        ));
    }

    // Assign the cost of graph nodes based on the input cost
    // vector. We also perform the "translation operation"
    // here to guaranttee a non-empty solution.
    transform_costs();

    // Build the arcs of the graphs.
    build_arcs();

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();

    // Recover optimal surface(s) from source set.
    for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
        for (i0 = 0; i0 < m_graph.size_0(); ++i0) {
            for (i2 = 0; i2 < m_graph.size_2(); ++i2) {
                
                // Find upper envelope.
                if (!m_graph.in_source_set(i0, i1, i2))
                    break;
            }
            net(i0, i1, 0) = (int)(i2 - 1);
        }
    }

    if (0 != pflow)
        *pflow = flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d<_Cost, _Cap, _Tg>::transform_costs()
{
    assert(m_pcost != 0);

    size_type            i0, i1, i2;

    const size_type& s0 = m_pcost->size_0();
    const size_type& s1 = m_pcost->size_1();
    const size_type& s2 = m_pcost->size_2();

    m_graph.set_initial_flow(0);

#   ifdef __OPTNET_SUPPORT_ROI__

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
    else { // else 1

#   endif // __OPTNET_SUPPORT_ROI__

        //
        // Construct the s-t graph "G_st".
        //
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

#   ifdef __OPTNET_SUPPORT_ROI__
    } // else 1
#   endif // __OPTNET_SUPPORT_ROI__

}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_fs_3d<_Cost, _Cap, _Tg>::build_arcs()
{
    int i0, i1, i2, s0, s1, s2;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();

    // Clear all constructed arcs.
    m_graph.clear_arcs();

    //
    // Construct graph arcs based on the given parameters.
    // -- Sorry, the lines here are a bit too long. ;-)
    //

    // Intra-column (vertical) arcs.
    for (i2 = s2 - 1; i2 > 0; --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                m_graph.add_arc(i0,          i1,          i2,          i0,          i1,          i2 - 1          );
            }
        }
    }

    // MKH -- only add arcs if more than one layer in this direction
    if (s0 > 1)
    {
    // Inter-column arcs (dir-0).
    for (i2 = s2 - 1; i2 > m_smooth[0]; --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 1; i0 < s0 - 1; ++i0) {
                m_graph.add_arc(i0,          i1,          i2,          i0 - 1,      i1,          i2 - m_smooth[0]);
                m_graph.add_arc(i0,          i1,          i2,          i0 + 1,      i1,          i2 - m_smooth[0]);
            } {
                m_graph.add_arc(0,           i1,          i2,          1,           i1,          i2 - m_smooth[0]);
                m_graph.add_arc(s0 - 1,      i1,          i2,          s0 - 2,      i1,          i2 - m_smooth[0]);
            }
        }
    }
    }
    /////

    // MKH -- only add arcs if more than one layer in this direction
    if (s1 > 1)
    {
    // Inter-column arcs (dir-1).
    for (i2 = s2 - 1; i2 > m_smooth[1]; --i2) {
        for (i0 = 0; i0 < s0; ++i0) {
            for (i1 = 1; i1 < s1 - 1; ++i1) {
                m_graph.add_arc(i0,          i1,          i2,          i0,          i1 - 1,      i2 - m_smooth[1]);
                m_graph.add_arc(i0,          i1,          i2,          i0,          i1 + 1,      i2 - m_smooth[1]);
            } {
                m_graph.add_arc(i0,          0,           i2,          i0,          1,           i2 - m_smooth[1]);
                m_graph.add_arc(i0,          s1 - 1,      i2,          i0,          s1 - 2,      i2 - m_smooth[1]);
            }
        }
    }
    }
    /////

    // The zero-plane.
    {
        // MKH -- only add arcs if more than one layer in this direction
        if (s0 > 1)
        {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 1; i0 < s0 - 1; ++i0) {
                m_graph.add_arc(i0,          i1,          0,           i0 - 1,      i1,          0               );
                m_graph.add_arc(i0,          i1,          0,           i0 + 1,      i1,          0               );
            } {
                m_graph.add_arc(0,           i1,          0,           1,           i1,          0               );
                m_graph.add_arc(s0 - 1,      i1,          0,           s0 - 2,      i1,          0               );
            }
        }
        }
        /////

        // MKH -- only add arcs if more than one layer in this direction
        if (s1 > 1)
        {
        for (i0 = 0; i0 < s0; ++i0) {
            for (i1 = 1; i1 < s1 - 1; ++i1) {
                m_graph.add_arc(i0,          i1,          0,           i0,          i1 - 1,      0               );
                m_graph.add_arc(i0,          i1,          0,           i0,          i1 + 1,      0               );
            } {
                m_graph.add_arc(i0,          0,           0,           i0,          1,           0               );
                m_graph.add_arc(i0,          s1 - 1,      0,           i0,          s1 - 2,      0               );
            }
        }
        }
        /////
    }

    // Circular graph connections.
    if (m_circle[0]) {
        for (i2 = s2 - 1; i2 > m_smooth[0]; --i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                m_graph.add_arc(0,           i1,          i2,          s0 - 1,      i1,          i2 - m_smooth[0]);
                m_graph.add_arc(s0 - 1,      i1,          i2,          0,           i1,          i2 - m_smooth[0]);
            }
        } { // i2 = 0
            for (i1 = 0; i1 < s1; ++i1) {
                m_graph.add_arc(0,           i1,          0,           s0 - 1,      i1,          0               );
                m_graph.add_arc(s0 - 1,      i1,          0,           0,           i1,          0               );
            }
        }
    } // if

    if (m_circle[1]) {
        for (i2 = s2 - 1; i2 > m_smooth[1]; --i2) {
            for (i0 = 0; i0 < s0; ++i0) {
                m_graph.add_arc(i0,          0,           i2,          i0,          s1 - 1,      i2 - m_smooth[1]);
                m_graph.add_arc(i0,          s1 - 1,      i2,          i0,          0,           i2 - m_smooth[1]);
            }
        } { // i2 = 0
            for (i0 = 0; i0 < s0; ++i0) {
                m_graph.add_arc(i0,          0,           0,           i0,          s1 - 1,      0               );
                m_graph.add_arc(i0,          s1 - 1,      0,           i0,          0,           0               );
            }
        }
    } // if
}

} // namespace

#endif
