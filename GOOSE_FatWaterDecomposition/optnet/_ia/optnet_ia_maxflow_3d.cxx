/*
 ==========================================================================
 |   
 |   $Id: optnet_ia_maxflow_3d.cxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___OPTNET_IA_MAXFLOW_3D_CXX___
#   define ___OPTNET_IA_MAXFLOW_3D_CXX___

#   include <optnet/_ia/optnet_ia_maxflow_3d.hxx>
#   include <limits>

#   ifdef max       // The max macro may interfere with
#       undef max   //   std::numeric_limits::max().
#   endif           //

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
optnet_ia_maxflow_3d<_Cap, _Tg>::optnet_ia_maxflow_3d() :
    m_proi(0), m_preflow(0)
{
    init();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
optnet_ia_maxflow_3d<_Cap, _Tg>::optnet_ia_maxflow_3d(size_type s0,
                                                      size_type s1,
                                                      size_type s2
                                                      ) :
    m_proi(0), m_preflow(0), m_nodes(s0, s1, s2)
{
    init();
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
bool
optnet_ia_maxflow_3d<_Cap, _Tg>::create(size_type s0,
                                        size_type s1,
                                        size_type s2
                                        )
{
    // Clears the ROI (Shouldn't do this!).
    // m_proi = 0;

    // Clear pre-calculated flow.
    m_preflow = 0;

    // Create and initialize node array.
    return m_nodes.create(s0, s1, s2);
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
typename optnet_ia_maxflow_3d<_Cap, _Tg>::capacity_type
optnet_ia_maxflow_3d<_Cap, _Tg>::solve()
{
    node_pointer   p_node, p_node1;
    node_pointer   p_node_cur     = 0;
    node_pointer   p_s_start_node = 0;
    node_pointer   p_t_start_node = 0;
    capacity_type* p_cap_rev_mid  = 0;
    capacity_type* p_cap_mid      = 0;
    int            i;
    
    // Initialize the maximum-flow solver.
    maxflow_init();

    while (true) {
        
        if (0 != (p_node = p_node_cur)) {
            p_node->tag &= ~IS_ACTIVE;
            if (!p_node->p_parent)
                p_node = 0;
        }

        if (!p_node) {

            while (!m_active_nodes.empty()) {
                p_node = m_active_nodes.front();
                p_node->tag &= ~IS_ACTIVE;
                m_active_nodes.pop_front();

                if (p_node->p_parent)
                    break;
            }
            if (m_active_nodes.empty())
                break;
        }
        
        //
        // Growth
        //
        p_s_start_node = 0;

        //
        // Growth -- Grow source tree.
        if (!(p_node->tag & IS_SINK)) {

            for (i = 0; i < m_max_arcs; ++i) {

                // Outgoing arcs.
                if (p_node->tag & m_mask_out[i]) {
                    p_node1 = p_node - m_offset[fwd_arc_index(p_node->tag)][i];
                }
                else p_node1 = 0;

                // Process adjacent nodes.
                if (0 != p_node1) {

                    if (!p_node1->p_parent) {
                        p_node1->tag       &= ~IS_SINK;
                        p_node1->tag       |= PARENT_REV;
                        p_node1->tag       &= ~PARENT_MASK;
                        p_node1->tag       |= (i & PARENT_MASK);
                        p_node1->p_parent   = p_node;
                        p_node1->dist_id    = p_node->dist_id;
                        p_node1->dist       = p_node->dist + 1;
                        activate(p_node1);
                    }
                    else if (p_node1->tag & IS_SINK) {
                        p_s_start_node      = p_node;
                        p_t_start_node      = p_node1;
                        p_cap_rev_mid       = &(p_node->res_cap[i]);
                        p_cap_mid           = 0;
                        break;
                    }
                    else if (p_node1->dist_id != 0 &&
                             p_node1->dist_id <= p_node->dist_id &&
                             p_node1->dist > p_node->dist) {
                        p_node1->tag       |= PARENT_REV;         //
                        p_node1->tag       &= ~PARENT_MASK;       //
                        p_node1->tag       |= (i & PARENT_MASK);  //
                        p_node1->p_parent   = p_node;             // Parent arc is reverse.
                        p_node1->dist_id    = p_node->dist_id;
                        p_node1->dist       = p_node->dist + 1;
                    }

                }
                
                // Incoming arcs
                if (p_node->tag & m_mask_in[i]) {
                    p_node1 = p_node + m_offset[rev_arc_index(p_node->tag)][i];
                }
                else p_node1 = 0;
                
                // Process adjacent nodes.
                if ((0 != p_node1) && (0 != p_node1->res_cap[i])) {

                    if (!p_node1->p_parent) {
                        p_node1->tag       &= ~IS_SINK;
                        p_node1->tag       &= ~PARENT_REV;        //
                        p_node1->tag       &= ~PARENT_MASK;       //
                        p_node1->tag       |= (i & PARENT_MASK);  //
                        p_node1->p_parent   = p_node;             // Parent arc is forward.
                        p_node1->dist_id    = p_node->dist_id;
                        p_node1->dist       = p_node->dist + 1;
                        activate(p_node1);
                    }
                    else if (p_node1->tag & IS_SINK) {
                        p_s_start_node      = p_node;
                        p_t_start_node      = p_node1;
                        p_cap_rev_mid        = 0;
                        p_cap_mid            = &(p_node1->res_cap[i]);
                        break;
                    }
                    else if (p_node1->dist_id != 0 &&
                             p_node1->dist_id <= p_node->dist_id && 
                             p_node1->dist > p_node->dist) {
                        p_node1->tag       &= ~PARENT_REV;        //
                        p_node1->tag       &= ~PARENT_MASK;       //
                        p_node1->tag       |= (i & PARENT_MASK);  //
                        p_node1->p_parent   = p_node;             // Parent arc is forward.
                        p_node1->dist_id    = p_node->dist_id;
                        p_node1->dist       = p_node->dist + 1;
                    }

                } // if (0 != p_node1) ...
            } // for (i = ...
        }
        //
        // Growth -- Grow sink tree.
        else {

            for (i = 0; i < m_max_arcs; ++i) {

                // Outgoing arcs
                if (p_node->tag & m_mask_out[i]) {
                    p_node1 = p_node - m_offset[fwd_arc_index(p_node->tag)][i];
                }
                else p_node1 = 0;
                
                // Process adjacent nodes.
                if ((0 != p_node1) && (0 != p_node->res_cap[i])) {

                    if (!p_node1->p_parent) {
                        p_node1->tag       |= IS_SINK;
                        p_node1->tag       |= PARENT_REV;         //
                        p_node1->tag       &= ~PARENT_MASK;       //
                        p_node1->tag       |= (i & PARENT_MASK);  //
                        p_node1->p_parent   = p_node;             // Parent arc is reverse.
                        p_node1->dist_id    = p_node->dist_id;
                        p_node1->dist       = p_node->dist + 1;
                        activate(p_node1);
                    }
                    else if (!(p_node1->tag & IS_SINK)) {
                        p_s_start_node      = p_node1;
                        p_t_start_node      = p_node;
                        p_cap_rev_mid        = 0;
                        p_cap_mid            = &(p_node->res_cap[i]);
                        break;
                    }
                    else if (p_node1->dist_id != 0 &&
                             p_node1->dist_id <= p_node->dist_id &&
                             p_node1->dist > p_node->dist) {
                        p_node1->tag       |= PARENT_REV;         //
                        p_node1->tag       &= ~PARENT_MASK;       //
                        p_node1->tag       |= (i & PARENT_MASK);  //
                        p_node1->p_parent   = p_node;             // Parent arc is reverse.
                        p_node1->dist_id    = p_node->dist_id;
                        p_node1->dist       = p_node->dist + 1;
                    }
                }

                // Incoming arcs.
                if (p_node->tag & m_mask_in[i]) {
                    p_node1 = p_node + m_offset[rev_arc_index(p_node->tag)][i];
                }
                else p_node1 = 0;
                
                // Process adjacent nodes.
                if (0 != p_node1) {
                    if (!p_node1->p_parent) {
                        p_node1->tag       |= IS_SINK;
                        p_node1->tag       &= ~PARENT_REV;        //
                        p_node1->tag       &= ~PARENT_MASK;       //
                        p_node1->tag       |= (i & PARENT_MASK);  //
                        p_node1->p_parent   = p_node;             // Parent arc is forward.
                        p_node1->dist_id    = p_node->dist_id;
                        p_node1->dist       = p_node->dist + 1;
                        activate(p_node1);
                    }
                    else if (!(p_node1->tag & IS_SINK)) {
                        p_s_start_node      = p_node1;
                        p_t_start_node      = p_node;
                        p_cap_rev_mid        = &(p_node1->res_cap[i]);
                        p_cap_mid            = 0;
                        break;
                    }
                    else if (p_node1->dist_id != 0 &&
                                p_node1->dist_id <= p_node->dist_id &&
                                p_node1->dist > p_node->dist) {
                        p_node1->tag       &= ~PARENT_REV;        //
                        p_node1->tag       &= ~PARENT_MASK;       //
                        p_node1->tag       |= (i & PARENT_MASK);  //
                        p_node1->p_parent   = p_node;             // Parent arc is forward.
                        p_node1->dist_id    = p_node->dist_id;
                        p_node1->dist       = p_node->dist + 1;
                    }
                }
            } // for (i = ...
        }

        if (0 != p_s_start_node) {

            p_node->tag |= IS_ACTIVE;
            p_node_cur   = p_node;
            
            maxflow_augment(p_s_start_node, 
                             p_t_start_node,
                             p_cap_mid,
                             p_cap_rev_mid);

            while (!m_orphan_nodes.empty()) {
                
                p_node = m_orphan_nodes.front();

                if (!(p_node->tag & IS_SINK))
                    maxflow_adopt_source_orphan(p_node);
                else
                    maxflow_adopt_sink_orphan(p_node);

                m_orphan_nodes.pop_front();
            }

            ++m_dist_id;
        }
        else
            p_node_cur = 0;

    } // while (true)

    return m_flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_ia_maxflow_3d<_Cap, _Tg>::init()
{
    m_mask_in [0] = 0x20000000L; // 0010 0000 0000 0000 ... -> x + 1
    m_mask_in [1] = 0x08000000L; // 0000 1000 0000 0000 ... -> x - 1
    m_mask_in [2] = 0x02000000L; // 0000 0010 0000 0000 ... -> y + 1
    m_mask_in [3] = 0x00800000L; // 0000 0000 1000 0000 ... -> y - 1
    m_mask_in [4] = 0x80000000L; // 1000 0000 0000 0000 ... -> z + 1

    m_mask_out[0] = 0x10000000L; // 0001 0000 0000 0000 ... -> x - 1
    m_mask_out[1] = 0x04000000L; // 0000 0100 0000 0000 ... -> x + 1
    m_mask_out[2] = 0x01000000L; // 0000 0001 0000 0000 ... -> y - 1
    m_mask_out[3] = 0x00400000L; // 0000 0000 0100 0000 ... -> y + 1
    m_mask_out[4] = 0x40000000L; // 0100 0000 0000 0000 ... -> z - 1

    // Default graph circlarty setting.
    m_circle[0] = false;
    m_circle[1] = false;

    // Default smoothness parameters.
    m_smooth[0] = 1;
    m_smooth[1] = 1;

    // Max num of arcs per node.
    m_max_arcs = 5;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_ia_maxflow_3d<_Cap, _Tg>::maxflow_init()
{
    int           i0, i1, i2, s0, s1, s2;
    node_pointer  p;

    s0 = (int)m_nodes.size_0();
    s1 = (int)m_nodes.size_1();
    s2 = (int)m_nodes.size_2();

    //
    // Initialize arc masks.
    //

    // Connect intra-column arcs.
    for (i2 = s2 - 1; i2 > 0; --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                if (in_roi (i0    , i1    , i2              ) &&
                    in_roi (i0    , i1    , i2 - 1          )) {
                    m_nodes(i0    , i1    , i2              ).tag |= m_mask_out[4];
                    m_nodes(i0    , i1    , i2 - 1          ).tag |= m_mask_in [4];
                }
            }
        }
    }

    // Connect inter-column arcs.
    for (i2 = s2 - 1; i2 > m_smooth[0]; --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 1; i0 < s0 - 1; ++i0) {
                if (in_roi (i0    , i1    , i2              ) && 
                    in_roi (i0 - 1, i1    , i2 - m_smooth[0])) {
                    m_nodes(i0    , i1    , i2              ).tag |= m_mask_out[0];
                    m_nodes(i0 - 1, i1    , i2 - m_smooth[0]).tag |= m_mask_in [0]; 
                }
                if (in_roi (i0    , i1    , i2              ) &&
                    in_roi (i0 + 1, i1    , i2 - m_smooth[0])) {
                    m_nodes(i0    , i1    , i2              ).tag |= m_mask_out[1];
                    m_nodes(i0 + 1, i1    , i2 - m_smooth[0]).tag |= m_mask_in [1]; 
                }
            } {
                if (in_roi (0     , i1    , i2              ) && 
                    in_roi (1     , i1    , i2 - m_smooth[0])) {
                    m_nodes(0     , i1    , i2              ).tag |= m_mask_out[1];
                    m_nodes(1     , i1    , i2 - m_smooth[0]).tag |= m_mask_in [1]; 
                }
                if (in_roi (s0 - 1, i1    , i2              ) && 
                    in_roi (s0 - 2, i1    , i2 - m_smooth[0])) {
                    m_nodes(s0 - 1, i1    , i2              ).tag |= m_mask_out[0]; 
                    m_nodes(s0 - 2, i1    , i2 - m_smooth[0]).tag |= m_mask_in [0]; 
                }
            }
            // Circular graph construction (x-direction).
            if (m_circle[0]) {
                // -- REMOVED --
                // if (i1 == 0 || i1 == (s1 - 1)) continue;
                if (in_roi (0     , i1    , i2              ) && 
                    in_roi (s0 - 1, i1    , i2 - m_smooth[0])) {
                    m_nodes(0     , i1    , i2              ).tag |= m_mask_out[0];
                    m_nodes(0     , i1    , i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(0     , i1    , i2              ).tag |= 0x080;
                    m_nodes(s0 - 1, i1    , i2 - m_smooth[0]).tag |= m_mask_in [0];
                    m_nodes(s0 - 1, i1    , i2 - m_smooth[0]).tag &= ~REV_ARC_MASK;
                    m_nodes(s0 - 1, i1    , i2 - m_smooth[0]).tag |= 0x800;
                }
                if (in_roi (s0 - 1, i1    , i2              ) && 
                    in_roi (0     , i1    , i2 - m_smooth[0])) {
                    m_nodes(s0 - 1, i1    , i2              ).tag |= m_mask_out[1];
                    m_nodes(s0 - 1, i1    , i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(s0 - 1, i1    , i2              ).tag |= 0x090;
                    m_nodes(0     , i1    , i2 - m_smooth[0]).tag |= m_mask_in [1];
                    m_nodes(0     , i1    , i2 - m_smooth[0]).tag &= ~REV_ARC_MASK;
                    m_nodes(0     , i1    , i2 - m_smooth[0]).tag |= 0x900;
                }
            } // if (m_circle[0])
        }
    }
    
    for (i2 = s2 - 1; i2 > m_smooth[1]; --i2) {
        for (i0 = 0; i0 < s0; ++i0) {
            for (i1 = 1; i1 < s1 - 1; ++i1) {
                if (in_roi (i0    , i1    , i2              ) &&
                    in_roi (i0    , i1 - 1, i2 - m_smooth[1])) {
                    m_nodes(i0    , i1    , i2              ).tag |= m_mask_out[2];
                    m_nodes(i0    , i1 - 1, i2 - m_smooth[1]).tag |= m_mask_in [2];
                }
                if (in_roi (i0    , i1    , i2              ) &&
                    in_roi (i0    , i1 + 1, i2 - m_smooth[1])) {
                    m_nodes(i0    , i1    , i2              ).tag |= m_mask_out[3];
                    m_nodes(i0    , i1 + 1, i2 - m_smooth[1]).tag |= m_mask_in [3];
                }
            } {                    
                if (in_roi (i0    , 0     , i2              ) && 
                    in_roi (i0    , 1     , i2 - m_smooth[1])) {
                    m_nodes(i0    , 0     , i2              ).tag |= m_mask_out[3];
                    m_nodes(i0    , 1     , i2 - m_smooth[1]).tag |= m_mask_in [3]; 
                }
                if (in_roi (i0    , s1 - 1, i2              ) &&
                    in_roi (i0    , s1 - 2, i2 - m_smooth[1])) {
                    m_nodes(i0    , s1 - 1, i2              ).tag |= m_mask_out[2]; 
                    m_nodes(i0    , s1 - 2, i2 - m_smooth[1]).tag |= m_mask_in [2]; 
                }
            }
            // Circular graph construction (y-direction).
            if (m_circle[1]) {
                // -- REMOVED --
                // if (i0 == 0 || i0 == (s0 - 1)) continue;
                if (in_roi (i0    , 0     , i2             ) &&
                    in_roi (i0    , s1 - 1, i2 - m_smooth[1])) {
                    m_nodes(i0    , 0     , i2              ).tag |= m_mask_out[2];
                    m_nodes(i0    , 0     , i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(i0    , 0     , i2              ).tag |= 0x060;
                    m_nodes(i0    , s1 - 1, i2 - m_smooth[1]).tag |= m_mask_in [2];
                    m_nodes(i0    , s1 - 1, i2 - m_smooth[1]).tag &= ~REV_ARC_MASK;
                    m_nodes(i0    , s1 - 1, i2 - m_smooth[1]).tag |= 0x600;
                }
                if (in_roi (i0    , s1 - 1, i2             ) &&
                    in_roi (i0    , 0     , i2 - m_smooth[1])) {
                    m_nodes(i0    , s1 - 1, i2              ).tag |= m_mask_out[3];
                    m_nodes(i0    , s1 - 1, i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(i0    , s1 - 1, i2              ).tag |= 0x070;
                    m_nodes(i0    , 0     , i2 - m_smooth[1]).tag |= m_mask_in [3];
                    m_nodes(i0    , 0     , i2 - m_smooth[1]).tag &= ~REV_ARC_MASK;
                    m_nodes(i0    , 0     , i2 - m_smooth[1]).tag |= 0x700;
                }
            } // if (m_circle[1])
        }
    }
    
    { // base set
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                {
                    m_nodes(i0    , i1    , lowest(i0    , i1   )).tag &= ~(FWD_ARC_MASK | REV_ARC_MASK);
                    m_nodes(i0    , i1    , lowest(i0    , i1   )).tag |= 0x110L;
                }
            }
        }

        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 1; i0 < s0 - 1; ++i0) {
                if (in_roi (i0    , i1    ) &&
                    in_roi (i0 - 1, i1    )) {
                    m_nodes(i0    , i1    , lowest(i0    , i1   )).tag |= m_mask_out[0];
                    m_nodes(i0 - 1, i1    , lowest(i0 - 1, i1   )).tag |= m_mask_in [0];
                }
                if (in_roi (i0    , i1    ) &&
                    in_roi (i0 + 1, i1    )) {
                    m_nodes(i0    , i1    , lowest(i0    , i1   )).tag |= m_mask_out[1];
                    m_nodes(i0 + 1, i1    , lowest(i0 + 1, i1   )).tag |= m_mask_in [1];
                }
            } {
                if (in_roi (0     , i1    ) &&
                    in_roi (1     , i1    )) {
                    m_nodes(0     , i1    , lowest(0     , i1   )).tag |= m_mask_out[1];
                    m_nodes(1     , i1    , lowest(1     , i1   )).tag |= m_mask_in [1]; 
                }
                if (in_roi (s0 - 1, i1    ) &&
                    in_roi (s0 - 2, i1    )) {
                    m_nodes(s0 - 1, i1    , lowest(s0 - 1, i1   )).tag |= m_mask_out[0]; 
                    m_nodes(s0 - 2, i1    , lowest(s0 - 2, i1   )).tag |= m_mask_in [0]; 
                }
            }
        }
        
        for (i0 = 0; i0 < s0; ++i0) {
            for (i1 = 1; i1 < s1 - 1; ++i1) {
                if (in_roi (i0    , i1    ) &&
                    in_roi (i0    , i1 - 1)) {
                    m_nodes(i0    , i1    , lowest(i0   , i1    )).tag |= m_mask_out[2];
                    m_nodes(i0    , i1 - 1, lowest(i0   , i1 - 1)).tag |= m_mask_in [2];
                }
                if (in_roi (i0    , i1    ) &&
                    in_roi (i0    , i1 + 1)) {
                    m_nodes(i0    , i1    , lowest(i0   , i1    )).tag |= m_mask_out[3];
                    m_nodes(i0    , i1 + 1, lowest(i0   , i1 + 1)).tag |= m_mask_in [3];
                }
            }
            if (s1 > 1) {                            
                if (in_roi (i0    , 0     ) &&
                    in_roi (i0    , 1     )) {
                    m_nodes(i0    , 0     , lowest(i0   , 0     )).tag |= m_mask_out[3];
                    m_nodes(i0    , 1     , lowest(i0   , 1     )).tag |= m_mask_in [3];
                }
                if (in_roi (i0    , s1 - 1) &&
                    in_roi (i0    , s1 - 2)) {
                    m_nodes(i0    , s1 - 1, lowest(i0   , s1 - 1)).tag |= m_mask_out[2];
                    m_nodes(i0    , s1 - 2, lowest(i0   , s1 - 2)).tag |= m_mask_in [2];
                }
            }
        }
    }

    if (m_circle[0] && m_circle[1]) {
        for (i2 = s2 - 1; i2 > m_smooth[0]; --i2) {
            if (s1 > 1) {
                if (in_roi (0     , 0     , i2              ) &&
                    in_roi (s0 - 1, 0     , i2 - m_smooth[0])) {
                    m_nodes(0     , 0     , i2              ).tag |= m_mask_out[0];
                    m_nodes(0     , 0     , i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(0     , 0     , i2              ).tag |= 0x020;
                    m_nodes(s0 - 1, 0     , i2 - m_smooth[0]).tag |= m_mask_in [0];
                    m_nodes(s0 - 1, 0     , i2 - m_smooth[0]).tag &= ~REV_ARC_MASK;
                    m_nodes(s0 - 1, 0     , i2 - m_smooth[0]).tag |= 0x400;//
                }
                if (in_roi (s0 - 1, 0     , i2              ) &&
                    in_roi (0     , 0     , i2 - m_smooth[0])) {
                    m_nodes(s0 - 1, 0     , i2              ).tag |= m_mask_out[1];
                    m_nodes(s0 - 1, 0     , i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(s0 - 1, 0     , i2              ).tag |= 0x050;
                    m_nodes(0     , 0     , i2 - m_smooth[0]).tag |= m_mask_in [1];
                    m_nodes(0     , 0     , i2 - m_smooth[0]).tag &= ~REV_ARC_MASK;
                    m_nodes(0     , 0     , i2 - m_smooth[0]).tag |= 0x300;//
                }

                if (in_roi (0     , s1 - 1, i2              ) &&
                    in_roi (s0 - 1, s1 - 1, i2 - m_smooth[0])) {
                    m_nodes(0     , s1 - 1, i2              ).tag |= m_mask_out[0];
                    m_nodes(0     , s1 - 1, i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(0     , s1 - 1, i2              ).tag |= 0x040;
                    m_nodes(s0 - 1, s1 - 1, i2 - m_smooth[0]).tag |= m_mask_in [0];
                    m_nodes(s0 - 1, s1 - 1, i2 - m_smooth[0]).tag &= ~REV_ARC_MASK;
                    m_nodes(s0 - 1, s1 - 1, i2 - m_smooth[0]).tag |= 0x200;//
                }
                if (in_roi (s0 - 1, s1 - 1, i2              ) &&
                    in_roi (0     , s1 - 1, i2 - m_smooth[0])) {
                    m_nodes(s0 - 1, s1 - 1, i2              ).tag |= m_mask_out[1];
                    m_nodes(s0 - 1, s1 - 1, i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(s0 - 1, s1 - 1, i2              ).tag |= 0x030;
                    m_nodes(0     , s1 - 1, i2 - m_smooth[0]).tag |= m_mask_in [1];
                    m_nodes(0     , s1 - 1, i2 - m_smooth[0]).tag &= ~REV_ARC_MASK;
                    m_nodes(0     , s1 - 1, i2 - m_smooth[0]).tag |= 0x500;//
                }
            } // if (s1 > 1)
        } // for (i2...
        for (i2 = s2 - 1; i2 > m_smooth[1]; --i2) {
            if (s0 > 1) {
                if (in_roi (0     , 0     , i2              ) &&
                    in_roi (0     , s1 - 1, i2 - m_smooth[1])) {
                    m_nodes(0     , 0     , i2              ).tag |= m_mask_out[2];
                    m_nodes(0     , 0     , i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(0     , 0     , i2              ).tag |= 0x020;
                    m_nodes(0     , s1 - 1, i2 - m_smooth[1]).tag |= m_mask_in [2];
                    m_nodes(0     , s1 - 1, i2 - m_smooth[1]).tag &= ~REV_ARC_MASK;
                    m_nodes(0     , s1 - 1, i2 - m_smooth[1]).tag |= 0x500;//
                }
                if (in_roi (0     , s1 - 1, i2              ) &&
                    in_roi (0     , 0     , i2 - m_smooth[1])) {
                    m_nodes(0     , s1 - 1, i2              ).tag |= m_mask_out[3];
                    m_nodes(0     , s1 - 1, i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(0     , s1 - 1, i2              ).tag |= 0x040;
                    m_nodes(0     , 0     , i2 - m_smooth[1]).tag |= m_mask_in [3];
                    m_nodes(0     , 0     , i2 - m_smooth[1]).tag &= ~REV_ARC_MASK;
                    m_nodes(0     , 0     , i2 - m_smooth[1]).tag |= 0x300;//
                }

                if (in_roi (s0 - 1, 0     , i2              ) &&
                    in_roi (s0 - 1, s1 - 1, i2 - m_smooth[1])) {
                    m_nodes(s0 - 1, 0     , i2              ).tag |= m_mask_out[2];
                    m_nodes(s0 - 1, 0     , i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(s0 - 1, 0     , i2              ).tag |= 0x050;
                    m_nodes(s0 - 1, s1 - 1, i2 - m_smooth[1]).tag |= m_mask_in [2];
                    m_nodes(s0 - 1, s1 - 1, i2 - m_smooth[1]).tag &= ~REV_ARC_MASK;
                    m_nodes(s0 - 1, s1 - 1, i2 - m_smooth[1]).tag |= 0x200;//
                }
                if (in_roi (s0 - 1, s1 - 1, i2              ) &&
                    in_roi (s0 - 1, 0     , i2 - m_smooth[1])) {
                    m_nodes(s0 - 1, s1 - 1, i2              ).tag |= m_mask_out[3];
                    m_nodes(s0 - 1, s1 - 1, i2              ).tag &= ~FWD_ARC_MASK;
                    m_nodes(s0 - 1, s1 - 1, i2              ).tag |= 0x030;
                    m_nodes(s0 - 1, 0     , i2 - m_smooth[1]).tag |= m_mask_in [3];
                    m_nodes(s0 - 1, 0     , i2 - m_smooth[1]).tag &= ~REV_ARC_MASK;
                    m_nodes(s0 - 1, 0     , i2 - m_smooth[1]).tag |= 0x400;//
                }
            } // if (s0 > 1)
        } // for (i2...
    }

    //
    // Initialize node offset table.
    //

    // These are for nodes not in the zero-plane.
    m_offset[0][0] = m_nodes.offset(1, 0, m_smooth[0]) - m_nodes.offset(0, 0, 0);       // x - 1, z - smooth0
    m_offset[0][1] = m_nodes.offset(0, 0, m_smooth[0]) - m_nodes.offset(1, 0, 0);       // x + 1, z - smooth0
    m_offset[0][2] = m_nodes.offset(0, 1, m_smooth[1]) - m_nodes.offset(0, 0, 0);       // y - 1, z - smooth1
    m_offset[0][3] = m_nodes.offset(0, 0, m_smooth[1]) - m_nodes.offset(0, 1, 0);       // y + 1, z - smooth1
    m_offset[0][4] = m_nodes.offset(0, 0, 1)           - m_nodes.offset(0, 0, 0);       // z - 1 

    // These are for nodes in the zero-plane.
    m_offset[1][0] = m_nodes.offset(1, 0, 0)           - m_nodes.offset(0, 0, 0);       // x - 1
    m_offset[1][1] = m_nodes.offset(0, 0, 0)           - m_nodes.offset(1, 0, 0);       // x + 1
    m_offset[1][2] = m_nodes.offset(0, 1, 0)           - m_nodes.offset(0, 0, 0);       // y - 1
    m_offset[1][3] = m_nodes.offset(0, 0, 0)           - m_nodes.offset(0, 1, 0);       // y + 1
    m_offset[1][4] = m_offset[0][4];                                                    // z - 1

    // These are for the circular graph.

    // These encode the different boundaries of the image.
    // See the figure below.
    //

    //
    // +--> x
    // |  
    // y  2o---6---+5
    //     |       |
    //     8       9
    //     |       |
    //    4+---7---+3
    //
    m_offset[2][0] = m_nodes.offset(0, 0, m_smooth[0]) - m_nodes.offset(s0 - 1, 0, 0);  // 0 -> s0-1
    m_offset[2][1] = m_offset[0][1];                                                    // x + 1, z - smooth0
    m_offset[2][2] = m_nodes.offset(0, 0, m_smooth[1]) - m_nodes.offset(0, s1 - 1, 0);  // 0 -> s1-1
    m_offset[2][3] = m_offset[0][3];                                                    // y + 1, z - smooth1
    m_offset[2][4] = m_offset[0][4];                                                    // z + 1

    m_offset[3][0] = m_offset[0][0];                                                    // x - 1, z - smooth0
    m_offset[3][1] = m_nodes.offset(s0 - 1, 0, m_smooth[0]) - m_nodes.offset(0, 0, 0);  // s0-1 -> 0
    m_offset[3][2] = m_offset[0][2];                                                    // y - 1, z - smooth1
    m_offset[3][3] = m_nodes.offset(0, s1 - 1, m_smooth[1]) - m_nodes.offset(0, 0, 0);  // s1-1 -> 0
    m_offset[3][4] = m_offset[0][4];                                                    // z + 1

    m_offset[4][0] = m_nodes.offset(0, 0, m_smooth[0]) - m_nodes.offset(s0 - 1, 0, 0);  // 0 -> s0-1
    m_offset[4][1] = m_offset[0][1];                                                    // x + 1, z - smooth0
    m_offset[4][2] = m_offset[0][2];                                                    // y - 1, z - smooth1
    m_offset[4][3] = m_nodes.offset(0, s1 - 1, m_smooth[1]) - m_nodes.offset(0, 0, 0);  // s1-1 -> 0
    m_offset[4][4] = m_offset[0][4];                                                    // z + 1

    m_offset[5][0] = m_offset[0][0];                                                    // x - 1, z - smooth0
    m_offset[5][1] = m_nodes.offset(s0 - 1, 0, m_smooth[0]) - m_nodes.offset(0, 0, 0);  // s0-1 -> 0
    m_offset[5][2] = m_nodes.offset(0, 0, m_smooth[1]) - m_nodes.offset(0, s1 - 1, 0);  // 0 -> s1-1
    m_offset[5][3] = m_offset[0][3];                                                    // y + 1, z - smooth1
    m_offset[5][4] = m_offset[0][4];                                                    // z + 1


    m_offset[6][0] = m_offset[0][0];                                                    // x - 1, z - smooth0
    m_offset[6][1] = m_offset[0][1];                                                    // x + 1, z - smooth0
    m_offset[6][2] = m_nodes.offset(0, 0, m_smooth[1]) - m_nodes.offset(0, s1 - 1, 0);  // 0 -> s1-1
    m_offset[6][3] = m_offset[0][3];                                                    // y + 1, z - smooth1
    m_offset[6][4] = m_offset[0][4];                                                    // z + 1

    m_offset[7][0] = m_offset[0][0];                                                    // x - 1, z - smooth0
    m_offset[7][1] = m_offset[0][1];                                                    // x + 1, z - smooth0
    m_offset[7][2] = m_offset[0][2];                                                    // y - 1, z - smooth1
    m_offset[7][3] = m_nodes.offset(0, s1 - 1, m_smooth[1]) - m_nodes.offset(0, 0, 0);  // s1-1 -> 0
    m_offset[7][4] = m_offset[0][4];                                                    // z + 1

    m_offset[8][0] = m_nodes.offset(0, 0, m_smooth[0]) - m_nodes.offset(s0 - 1, 0, 0);  // 0 -> s0-1
    m_offset[8][1] = m_offset[0][1];                                                    // x + 1, z - smooth0
    m_offset[8][2] = m_offset[0][2];                                                    // y - 1, z - smooth1
    m_offset[8][3] = m_offset[0][3];                                                    // y + 1, z - smooth1
    m_offset[8][4] = m_offset[0][4];                                                    // z + 1

    m_offset[9][0] = m_offset[0][0];                                                    // x - 1, z - smooth0
    m_offset[9][1] = m_nodes.offset(s0 - 1, 0, m_smooth[0]) - m_nodes.offset(0, 0, 0);  // s0-1 -> 0
    m_offset[9][2] = m_offset[0][2];                                                    // y - 1, z - smooth1
    m_offset[9][3] = m_offset[0][3];                                                    // y + 1, z - smooth1
    m_offset[9][4] = m_offset[0][4];                                                    // z + 1

    //
    // Initialize node queues and flow.
    //
    m_active_nodes.clear();
    m_orphan_nodes.clear();
    m_flow = m_preflow;
    
    for (p = &*m_nodes.begin(); p != &*m_nodes.end(); ++p) {

        if (p->cap > 0) {
            // The node is connected to the source.
            p->tag     &= ~IS_SINK;
            p->p_parent = (node_pointer)TERMINAL;
            p->dist_id  = 1;
            p->dist     = 1;
            activate(p);
        }
        else if (p->cap < 0) {
            // The node is connected to the sink
            p->tag     |= IS_SINK;
            p->p_parent = (node_pointer)TERMINAL;
            p->dist_id  = 1;
            p->dist     = 1;
            activate(p);
        }
        else {
            p->p_parent = 0;
            p->dist_id  = 0;
            p->dist    = 0;
        }
    }

    //
    // Global distance ID.
    //
    m_dist_id = 2;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_ia_maxflow_3d<_Cap, _Tg>::maxflow_augment(
                                    node_pointer   p_s_start_node,
                                    node_pointer   p_t_start_node, 
                                    capacity_type* p_cap_mid, 
                                    capacity_type* p_cap_rev_mid
                                    )
{

    capacity_type   bottle_neck_cap;
    node_pointer    p_node, p_node1;
    int             i;

    // STEP 1: find bottleneck capacity
    bottle_neck_cap = p_cap_mid ? 
        *p_cap_mid : std::numeric_limits<capacity_type>::max();

    //  1-1: the source tree
    for (p_node = p_s_start_node; 
         p_node->p_parent != (node*)TERMINAL;
         p_node = p_node->p_parent) {

        if (!(p_node->tag & PARENT_REV)) {
            i = (int)(p_node->tag & PARENT_MASK);
            if (bottle_neck_cap > p_node->res_cap[i])
                bottle_neck_cap = p_node->res_cap[i];
        }
    }

    if (bottle_neck_cap > p_node->cap)
        bottle_neck_cap = p_node->cap;
    if (bottle_neck_cap == 0) return;

    //  1-2: the sink tree
    for (p_node = p_t_start_node;
         p_node->p_parent != (node*)TERMINAL;
         p_node = p_node->p_parent) {

        if (p_node->tag & PARENT_REV) {

            i = (int)(p_node->tag & PARENT_MASK);

            if (bottle_neck_cap > p_node->p_parent->res_cap[i])
                bottle_neck_cap = p_node->p_parent->res_cap[i];
        }
    }
    
    if (bottle_neck_cap > -p_node->cap)
        bottle_neck_cap = -p_node->cap;
    if (bottle_neck_cap == 0) return;


    // STEP 2: augment
    if (p_cap_rev_mid) *p_cap_rev_mid += bottle_neck_cap;
    if (p_cap_mid)     *p_cap_mid     -= bottle_neck_cap;

    //  2-1: the source tree
    for (p_node = p_s_start_node;
         p_node->p_parent != (node*)TERMINAL;
         p_node = p_node1) {

        p_node1 = p_node->p_parent;

        if (!(p_node->tag & PARENT_REV)) {

            i = (int)(p_node->tag & PARENT_MASK);
            
            p_node->res_cap[i] -= bottle_neck_cap;
            
            if (!p_node->res_cap[i]) {
                p_node->p_parent = (node*)ORPHAN;
                m_orphan_nodes.push_back(p_node);
            }
        }
        else {
            i = (int)(p_node->tag & PARENT_MASK);
            p_node1->res_cap[i] += bottle_neck_cap;
        } // else

    } // for

    p_node->cap -= bottle_neck_cap;
    if (!p_node->cap) {
        p_node->p_parent = (node*)ORPHAN;
        m_orphan_nodes.push_back(p_node);
    }

    //  2-2: the sink tree
    for (p_node = p_t_start_node;
         p_node->p_parent != (node*)TERMINAL;
         p_node = p_node1) {

        p_node1 = p_node->p_parent;

        if (p_node->tag & PARENT_REV) {

            i = (int)(p_node->tag & PARENT_MASK);

            p_node1->res_cap[i] -= bottle_neck_cap;

            if (!p_node1->res_cap[i]) {
                p_node->p_parent = (node*)ORPHAN;
                m_orphan_nodes.push_back(p_node);
            }
        }
        else {
            i = (int)(p_node->tag & PARENT_MASK);
            p_node->res_cap[i] += bottle_neck_cap;
        } // else

    } // for

    p_node->cap += bottle_neck_cap;
    if (!p_node->cap) {
        p_node->p_parent = (node*)ORPHAN;
        m_orphan_nodes.push_back(p_node);
    }

    m_flow += bottle_neck_cap;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_ia_maxflow_3d<_Cap, _Tg>::maxflow_adopt_source_orphan(
                                    node_pointer p_orphan
                                    )
{
    static const size_type
        DIST_MAX = std::numeric_limits<size_type>::max();

    struct parent_info_ {
        node_pointer  p_node;
        int           arc;
        bool          is_rev;
    } min_parent   = { 0, -1, false };
    node_pointer     p_node, p_node1, p_node2;
    size_type        dist, min_dist = DIST_MAX;
    int              i;


    // Try to find new parents for the orphan.
    for (i = 0; i < m_max_arcs; ++i) {
        if (0 != p_orphan->res_cap[i]) { // if 1
            //
            // Outgoing arcs.
            if (p_orphan->tag & m_mask_out[i]) {
                p_node = p_orphan - m_offset[fwd_arc_index(p_orphan->tag)][i];
            }
            else p_node = 0;
            
            if ((0 != p_node) && 
                (!(p_node->tag & IS_SINK)) && 
                (p_node->p_parent)) { // if 2

                p_node1 = p_node;   // Save node iterator.
                dist     = 0;         // Reset distance.

                while (true) {

                    if (p_node1->dist_id == m_dist_id) {
                        dist += p_node1->dist;
                        break;
                    }
                    ++dist;

                    p_node2 = p_node1->p_parent;

                    if (p_node2 == (node*)TERMINAL) {
                        p_node1->dist_id = m_dist_id;
                        p_node1->dist = 1;
                        break;
                    } 
                    
                    if (p_node2 == (node*)ORPHAN) {
                        dist = DIST_MAX;
                        break;
                    }

                    // Trace back to the source.
                    p_node1 = p_node2;

                } // while

                if (dist < DIST_MAX) {

                    // Save minimum distance node so far.
                    if (dist < min_dist) {
                        min_parent.p_node = p_node;
                        min_parent.arc    = i;
                        min_parent.is_rev  = false;
                        min_dist           = dist;
                    }
                    
                    // Set distance along the path.
                    for (;
                        p_node->dist_id != m_dist_id;
                        p_node = p_node->p_parent) {

                        p_node->dist_id = m_dist_id;
                        p_node->dist = dist--;
                    }
                }
            } // if 2
        } // if 1

        //
        // Incoming arcs.
        if (p_orphan->tag & m_mask_in[i]) {
            p_node = p_orphan + m_offset[rev_arc_index(p_orphan->tag)][i];
        }
        else p_node = 0;

        if ((0 != p_node) && 
            (!(p_node->tag & IS_SINK)) && 
            (p_node->p_parent)) { // if 3

            p_node1 = p_node;   // Save node iterator.
            dist     = 0;         // Reset distance.

            while (true) {

                if (p_node1->dist_id == m_dist_id) {
                    dist += p_node1->dist;
                    break;
                }
                ++dist;

                p_node2 = p_node1->p_parent;

                if (p_node2 == (node*)TERMINAL) {
                    p_node1->dist_id = m_dist_id;
                    p_node1->dist = 1;
                    break;
                } 
                
                if (p_node2 == (node*)ORPHAN) {
                    dist = DIST_MAX;
                    break;
                }

                // Trace back to the source.
                p_node1 = p_node2;

            } // while

            if (dist < DIST_MAX) {

                // Save minimum distance node so far.
                if (dist < min_dist) {
                    min_parent.p_node = p_node;
                    min_parent.arc    = i;
                    min_parent.is_rev  = true;
                    min_dist           = dist;
                }
                
                // Set distance along the path.
                for (;
                    p_node->dist_id != m_dist_id;
                    p_node = p_node->p_parent) {

                    p_node->dist_id = m_dist_id;
                    p_node->dist = dist--;
                }
            }
        } // if 3

    } // for

    // If found parent.
    if (0 != (p_orphan->p_parent = min_parent.p_node)) {

        if (min_parent.is_rev) p_orphan->tag |=  PARENT_REV;
        else                   p_orphan->tag &= ~PARENT_REV;

        p_orphan->tag &= ~PARENT_MASK;
        p_orphan->tag |= (PARENT_MASK & min_parent.arc);

        p_orphan->dist_id = m_dist_id;
        p_orphan->dist = min_dist + 1;

    } else { // Parent not found.
        
        p_orphan->dist_id = 0;

        // Process neighbors.
        for (i = 0; i < m_max_arcs; ++i) {
            //
            // Outgoing arcs.
            if (p_orphan->tag & m_mask_out[i]) {
                p_node = p_orphan - m_offset[fwd_arc_index(p_orphan->tag)][i];
            }
            else 
                p_node = 0;
            
            if ((0 != p_node) && 
                (!(p_node->tag & IS_SINK)) && 
                (0 != (p_node1 = p_node->p_parent))) { // if 4
                
                if (0 != p_orphan->res_cap[i]) activate(p_node);

                if (p_node1 == p_orphan) {
                    p_node->p_parent = (node*)ORPHAN;
                    m_orphan_nodes.push_back(p_node);
                }
            } // if

            //
            // Incoming arcs.
            if (p_orphan->tag & m_mask_in[i]) {
                p_node = p_orphan + m_offset[rev_arc_index(p_orphan->tag)][i];
            }
            else
                p_node = 0;

            if ((0 != p_node) && 
                (!(p_node->tag & IS_SINK)) && 
                (0 != (p_node1 = p_node->p_parent))) { // if 4

                activate(p_node);

                if (p_node1 == p_orphan) {
                    p_node->p_parent = (node*)ORPHAN;
                    m_orphan_nodes.push_back(p_node);
                }
            } // if
        } // for
    } // else
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_ia_maxflow_3d<_Cap, _Tg>::maxflow_adopt_sink_orphan(
                                    node_pointer p_orphan
                                    )
{
    static const size_type
        DIST_MAX = std::numeric_limits<size_type>::max();

    struct parent_info_ {
        node_pointer  p_node;
        int           arc;
        bool          is_rev;
    } min_parent   = { 0, -1, false };
    node_pointer     p_node, p_node1, p_node2;
    size_type        dist, min_dist = DIST_MAX;
    int              i;

    // Try to find new parents for the orphan.
    for (i = 0; i < m_max_arcs; ++i) {

        //
        // Outgoing arcs.
        if (p_orphan->tag & m_mask_out[i]) {
            p_node = p_orphan - m_offset[fwd_arc_index(p_orphan->tag)][i];
        }
        else
            p_node = 0;
        
        if ((0 != p_node) && 
            (p_node->tag & IS_SINK) && 
            (p_node->p_parent)) { // if 1

            p_node1 = p_node;   // Save node iterator.
            dist     = 0;         // Reset distance.

            while (true) {

                if (p_node1->dist_id == m_dist_id) {
                    dist += p_node1->dist;
                    break;
                }
                ++dist;

                p_node2 = p_node1->p_parent;

                if (p_node2 == (node*)TERMINAL) {
                    p_node1->dist_id = m_dist_id;
                    p_node1->dist = 1;
                    break;
                } 
                
                if (p_node2 == (node*)ORPHAN) {
                    dist = DIST_MAX;
                    break;
                }

                // Trace back to the source.
                p_node1 = p_node2;

            } // while

            if (dist < DIST_MAX) {

                // Save minimum distance node so far.
                if (dist < min_dist) {
                    min_parent.p_node = p_node;
                    min_parent.arc     = i;
                    min_parent.is_rev  = false;
                    min_dist           = dist;
                }
                
                // Set distance along the path.
                for (;
                    p_node->dist_id != m_dist_id;
                    p_node = p_node->p_parent) {

                    p_node->dist_id = m_dist_id;
                    p_node->dist = dist--;
                }
            }

        } // if 1
        
        //
        // Incoming arcs
        if (p_orphan->tag & m_mask_in[i]) {
            p_node = p_orphan + m_offset[rev_arc_index(p_orphan->tag)][i];
        }
        else
            p_node = 0;
            
        if ((0 != p_node) && 
            (0 != p_node->res_cap[i]) && 
            (p_node->tag & IS_SINK) && 
            (p_node->p_parent)) { // if 2

            p_node1 = p_node;   // Save node iterator.
            dist     = 0;         // Reset distance.

            while (true) {

                if (p_node1->dist_id == m_dist_id) {
                    dist += p_node1->dist;
                    break;
                }
                ++dist;

                p_node2 = p_node1->p_parent;

                if (p_node2 == (node*)TERMINAL) {
                    p_node1->dist_id = m_dist_id;
                    p_node1->dist = 1;
                    break;
                } 
                
                if (p_node2 == (node*)ORPHAN) {
                    dist = DIST_MAX;
                    break;
                }

                // Trace back to the source.
                p_node1 = p_node2;

            } // while

            if (dist < DIST_MAX) {

                // Save minimum distance node so far.
                if (dist < min_dist) {
                    min_parent.p_node = p_node;
                    min_parent.arc     = i;
                    min_parent.is_rev  = true;
                    min_dist           = dist;
                }
                
                // Set distance along the path.
                for (;
                    p_node->dist_id != m_dist_id;
                    p_node = p_node->p_parent) {

                    p_node->dist_id = m_dist_id;
                    p_node->dist = dist--;
                }
            }
        } // if 2
    } // for

    // If found parent.
    if (0 != (p_orphan->p_parent = min_parent.p_node)) {
        
        if (min_parent.is_rev) p_orphan->tag |=  PARENT_REV;
        else                   p_orphan->tag &= ~PARENT_REV;

        p_orphan->tag &= ~PARENT_MASK;
        p_orphan->tag |= (PARENT_MASK & min_parent.arc);

        p_orphan->dist_id = m_dist_id;
        p_orphan->dist = min_dist + 1;

    }
    // Parent not found.
    else {

        p_orphan->dist_id = 0;

        // Process neighbors.
        for (i = 0; i < m_max_arcs; ++i) {

            //
            // Outgoing arcs.
            if (p_orphan->tag & m_mask_out[i]) {
                p_node = p_orphan - m_offset[fwd_arc_index(p_orphan->tag)][i];
            }
            else
                p_node = 0;
            
            if ((0 != p_node) && 
                (p_node->tag & IS_SINK) && 
                (0 != (p_node1 = p_node->p_parent))) { // if 3

                activate(p_node);

                if (p_node1 == p_orphan) {
                    p_node->p_parent = (node*)ORPHAN;
                    m_orphan_nodes.push_back(p_node);
                }
            } // if 3

            //
            // Incoming arcs.
            if (p_orphan->tag & m_mask_in[i]) {
                p_node = p_orphan + m_offset[rev_arc_index(p_orphan->tag)][i];
            }
            else
                p_node = 0;

            if ((0 != p_node) && 
                (p_node->tag & IS_SINK) && 
                (0 != (p_node1 = p_node->p_parent))) { // if 4

                if (0 != p_node->res_cap[i]) activate(p_node);

                if (p_node1 == p_orphan) {
                    p_node->p_parent = (node*)ORPHAN;
                    m_orphan_nodes.push_back(p_node);
                }
            } // if 4
        } // for
    } // else
}

//
// Constants
//
template <typename _Cap, typename _Tg>
    const typename optnet_ia_maxflow_3d<_Cap, _Tg>::node* 
        optnet_ia_maxflow_3d<_Cap, _Tg>::TERMINAL =
            (typename optnet_ia_maxflow_3d<_Cap, _Tg>::node*)(1);

template <typename _Cap, typename _Tg>
    const typename optnet_ia_maxflow_3d<_Cap, _Tg>::node* 
        optnet_ia_maxflow_3d<_Cap, _Tg>::ORPHAN =
            (typename optnet_ia_maxflow_3d<_Cap, _Tg>::node*)(2);

template <typename _Cap, typename _Tg>
    const long optnet_ia_maxflow_3d<_Cap, _Tg>::FWD_ARC_MASK = 0x000000F0L;

template <typename _Cap, typename _Tg>
    const long optnet_ia_maxflow_3d<_Cap, _Tg>::REV_ARC_MASK = 0x00000F00L;

template <typename _Cap, typename _Tg>
    const long optnet_ia_maxflow_3d<_Cap, _Tg>::PARENT_MASK  = 0x0000000FL;

template <typename _Cap, typename _Tg>
    const long optnet_ia_maxflow_3d<_Cap, _Tg>::PARENT_REV   = 0x00080000L;

template <typename _Cap, typename _Tg>
    const long optnet_ia_maxflow_3d<_Cap, _Tg>::IS_ACTIVE    = 0x00020000L;

template <typename _Cap, typename _Tg>
    const long optnet_ia_maxflow_3d<_Cap, _Tg>::IS_SINK      = 0x00010000L;

} // namespace

#endif
