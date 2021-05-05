/*
 ==========================================================================
 |   
 |   $Id: optnet_fs_maxflow.cxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___OPTNET_FS_MAXFLOW_CXX___
#   define ___OPTNET_FS_MAXFLOW_CXX___

#   include <optnet/_fs/optnet_fs_maxflow.hxx>
#   include <limits>

#   ifdef max       // The max macro may interfere with
#       undef max   //   std::numeric_limits::max().
#   endif           //

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
optnet_fs_maxflow<_Cap, _Tg>::optnet_fs_maxflow() :
    m_preflow(0)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
optnet_fs_maxflow<_Cap, _Tg>::optnet_fs_maxflow(size_type s0,
                                                size_type s1,
                                                size_type s2,
                                                size_type s3,
                                                size_type s4
                                                ) :
    _Base(s0, s1, s2, s3, s4), m_preflow(0)
{

}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
bool
optnet_fs_maxflow<_Cap, _Tg>::create(size_type s0,
                                     size_type s1,
                                     size_type s2,
                                     size_type s3,
                                     size_type s4
                                     )
{
    // Clear pre-calculated flow.
    m_preflow = 0;

    // Create and initialize the graph.
    return _Base::create(s0, s1, s2, s3, s4);
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
typename optnet_fs_maxflow<_Cap, _Tg>::capacity_type
optnet_fs_maxflow<_Cap, _Tg>::solve()
{
    forward_arc_pointer p_fwd_arc, p_first_out_arc, p_last_out_arc;
    reverse_arc_pointer p_rev_arc, p_first_in_arc, p_last_in_arc;
    node_pointer        p_s_start_node = 0, p_t_start_node = 0;
    node_pointer        p_node, p_node1, p_cur_node = 0;
    capacity_type*      p_mid_rev_cap = 0;
    capacity_type*      p_mid_fwd_cap = 0;
    

    // Initialize the maximum-flow solver.
    maxflow_init();

    while (true) {

        if (0 != (p_node = p_cur_node)) {

            p_node->tag &= ~IS_ACTIVE;
            if (!p_node->p_parent_arc)
                p_node = 0;
        }
        if (!p_node) {

            while (!m_active_nodes.empty()) {

                p_node = m_active_nodes.front();
                p_node->tag &= ~IS_ACTIVE;
                m_active_nodes.pop_front();

                if (p_node->p_parent_arc)
                    break;
            }

            if (m_active_nodes.empty())
                break;
        }

        //
        // STAGE 1: Growth.
        //
        p_s_start_node = 0;
        
        ///////////////////////////////////////////////////////////////////
        // Normally, p_node->p_first_out_arc is the first outgoing arc of
        // p_node and (p_node+1)->p_first_out_arc-1 is the last outgoing
        // arc.
        //
        // However, since the arcs are stored in chunks, the arcs incident
        // to two consecutive nodes may be in different chunks.
        //
        // To handle this,  different schemes are used for the nodes whose
        // arcs are last in an arc chunk. In this case, the first outgoing
        // arc of p_node is p_node->p_first_out_arc+1, and the last out-
        // going arc is p_node->p_first_out_arc->shift.  Accordingly, the
        // first incoming arc of p_node is p_node->p_first_in_arc+1, and
        // the last incoming arc is p_node->p_first_in_arc->p_fwd.
        //
        if (p_node->tag & _Base::SPECIAL_OUT) {
            p_first_out_arc = p_node->p_first_out_arc + 1;
            p_last_out_arc
                = (forward_arc_pointer)(p_node->p_first_out_arc->shift);
        }
        else {
            p_first_out_arc = p_node->p_first_out_arc;
            p_last_out_arc  = (p_node + 1)->p_first_out_arc;
        }
        
        if (p_node->tag & _Base::SPECIAL_IN) {
            p_first_in_arc  = p_node->p_first_in_arc + 1;
            p_last_in_arc
                = (reverse_arc_pointer)(p_node->p_first_in_arc->p_fwd);
        }
        else {
            p_first_in_arc  = p_node->p_first_in_arc;
            p_last_in_arc   = (p_node + 1)->p_first_in_arc;
        }
        
        // Grow source tree.
        if (!(p_node->tag & IS_SINK)) {

            for (p_fwd_arc = p_first_out_arc;
                 p_fwd_arc < p_last_out_arc;
                 ++p_fwd_arc) {

                p_node1 = neighbor_node_fwd(p_node, p_fwd_arc->shift);

                if (!p_node1->p_parent_arc) {
                    // Parent arc is reverse.
                    p_node1->tag          &= ~IS_SINK;
                    p_node1->tag          |= PARENT_REV;
                    p_node1->p_parent_arc  = p_fwd_arc;
                    p_node1->dist_id       = p_node->dist_id;
                    p_node1->dist          = p_node->dist + 1;
                    activate(p_node1);
                }
                else if (p_node1->tag & IS_SINK) {
                    p_s_start_node         = p_node;
                    p_t_start_node         = p_node1;
                    p_mid_fwd_cap          = 0;
                    p_mid_rev_cap          = &(p_fwd_arc->rev_cap);
                    break;
                }
                else if (p_node1->dist_id != 0 &&
                         p_node1->dist_id <= p_node->dist_id &&
                         p_node1->dist > p_node->dist) {
                    // Parent arc is reverse.
                    p_node1->tag          |= PARENT_REV;
                    p_node1->p_parent_arc  = p_fwd_arc;
                    p_node1->dist_id       = p_node->dist_id;
                    p_node1->dist          = p_node->dist + 1;
                }

            } // for (p_fwd_arc ...

            if (!p_s_start_node) {

                for (p_rev_arc = p_first_in_arc;
                     p_rev_arc < p_last_in_arc;
                     ++p_rev_arc) {

                    p_fwd_arc = p_rev_arc->p_fwd;

                    if (0 != p_fwd_arc->rev_cap) {

                        p_node1 = neighbor_node_rev(p_node, p_fwd_arc->shift);

                        if (!p_node1->p_parent_arc) {
                            // Parent arc is forward.
                            p_node1->tag          &= ~IS_SINK;
                            p_node1->tag          &= ~PARENT_REV;
                            p_node1->p_parent_arc  = p_fwd_arc;
                            p_node1->dist_id       = p_node->dist_id;
                            p_node1->dist          = p_node->dist + 1;
                            activate(p_node1);
                        }
                        else if (p_node1->tag & IS_SINK) {
                            p_s_start_node         = p_node;
                            p_t_start_node         = p_node1;
                            p_mid_fwd_cap           = &(p_fwd_arc->rev_cap);
                            p_mid_rev_cap           = 0;
                            break;
                        }
                        else if (p_node1->dist_id != 0 &&
                                 p_node1->dist_id <= p_node->dist_id &&
                                 p_node1->dist > p_node->dist) {
                            // Parent arc is forward.
                            p_node1->tag          &= ~PARENT_REV;
                            p_node1->p_parent_arc  = p_fwd_arc;
                            p_node1->dist_id       = p_node->dist_id;
                            p_node1->dist          = p_node->dist + 1;
                        }
                    } // if (0 != p_fwd_arc->rev_cap)
                } // for (p_rev_arc ...
            } // if (!p_s_start_node)
        }
        else {
            // Grow sink tree.
            for (p_fwd_arc = p_first_out_arc;
                 p_fwd_arc < p_last_out_arc;
                 ++p_fwd_arc) {

                if (0 != p_fwd_arc->rev_cap) {

                    p_node1 = neighbor_node_fwd(p_node, p_fwd_arc->shift);

                    if (!p_node1->p_parent_arc) {
                        // Parent arc is reverse.
                        p_node1->tag          |= IS_SINK;
                        p_node1->tag          |= PARENT_REV;
                        p_node1->p_parent_arc  = p_fwd_arc;        
                        p_node1->dist_id       = p_node->dist_id;
                        p_node1->dist          = p_node->dist + 1;
                        activate(p_node1);
                    }
                    else if (!(p_node1->tag & IS_SINK)) {
                        p_s_start_node         = p_node1;
                        p_t_start_node         = p_node;
                        p_mid_fwd_cap           = &(p_fwd_arc->rev_cap);
                        p_mid_rev_cap           = 0;
                        break;
                    }
                    else if (p_node1->dist_id != 0 &&
                             p_node1->dist_id <= p_node->dist_id &&
                             p_node1->dist > p_node->dist) {
                        // Parent arc is reverse.
                        p_node1->tag          |= PARENT_REV;
                        p_node1->p_parent_arc  = p_fwd_arc;        
                        p_node1->dist_id       = p_node->dist_id;
                        p_node1->dist          = p_node->dist + 1;
                    }
                } // if (0 != p_fwd_arc->rev_cap)
            } // for (p_fwd_arc ...

            if (!p_s_start_node) {

                for (p_rev_arc = p_first_in_arc;
                     p_rev_arc < p_last_in_arc;
                     ++p_rev_arc) {

                    p_fwd_arc = p_rev_arc->p_fwd;

                    p_node1 = neighbor_node_rev(p_node, p_fwd_arc->shift);

                    if (!p_node1->p_parent_arc) {
                        // Parent arc is forward.
                        p_node1->tag          |= IS_SINK;
                        p_node1->tag          &= ~PARENT_REV;      
                        p_node1->p_parent_arc  = p_fwd_arc;        
                        p_node1->dist_id       = p_node->dist_id;
                        p_node1->dist          = p_node->dist + 1;
                        activate(p_node1);
                    }
                    else if (!(p_node1->tag & IS_SINK)) {
                        p_s_start_node         = p_node1;
                        p_t_start_node         = p_node;
                        p_mid_fwd_cap           = 0;
                        p_mid_rev_cap           = &(p_fwd_arc->rev_cap);
                        break;
                    }
                    else if (p_node1->dist_id != 0 &&
                             p_node1->dist_id <= p_node->dist_id &&
                             p_node1->dist > p_node->dist) {
                        // Parent arc is forward.
                        p_node1->tag          &= ~PARENT_REV;
                        p_node1->p_parent_arc  = p_fwd_arc;        
                        p_node1->dist_id       = p_node->dist_id;
                        p_node1->dist          = p_node->dist + 1;
                    }
                
                } // for (p_rev_arc ...
            } // if (!p_s_start_node)
        } // else

        if (0 == p_s_start_node) p_cur_node = 0;
        else {

            // Set active flag (prevent reactivation).
            p_node->tag |= IS_ACTIVE;
            p_cur_node   = p_node;
            
            //
            // STAGE 2: Augmentation.
            //
            maxflow_augment(p_s_start_node,
                            p_t_start_node,
                            p_mid_fwd_cap,
                            p_mid_rev_cap
                            );

            //
            // STAGE 3: Adoption.
            //
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

    } // while (true)

    return m_flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_fs_maxflow<_Cap, _Tg>::maxflow_init()
{
    node_pointer p;

    // Convert the graph to a forward-star representation.
    _Base::prepare();

    // Initialize data structures.
    m_active_nodes.clear();
    m_orphan_nodes.clear();
    m_flow = m_preflow;

    for (p = &*_Base::m_nodes.begin(); p != &*_Base::m_nodes.end(); ++p) {

        if (p->cap > 0) {

            // p is connected to the source
            p->tag          &= ~IS_SINK;
            p->p_parent_arc  = (forward_arc_pointer)TERMINAL;
            p->dist_id       = 1;
            p->dist          = 1;
            activate(p);
        }
        else if (p->cap < 0) {

            // p is connected to the sink
            p->tag          |= IS_SINK;
            p->p_parent_arc  = (forward_arc_pointer)TERMINAL;
            p->dist_id       = 1;
            p->dist          = 1;
            activate(p);
        }
        else {
            p->p_parent_arc  = NULL;
            p->dist_id       = 0;
            p->dist          = 0;
        }

    } // for

    m_dist_id = 2;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_fs_maxflow<_Cap, _Tg>::maxflow_augment(
                                node_pointer   p_s_start_node,
                                node_pointer   p_t_start_node,
                                capacity_type* p_mid_fwd_cap,
                                capacity_type* p_mid_rev_cap
                                )
{
    capacity_type       bottle_neck_cap;
    forward_arc_pointer p_fwd_arc;
    node_pointer        p_node;

    // STEP 1: Find bottleneck capacity.
    bottle_neck_cap 
        = p_mid_fwd_cap ? 
            *p_mid_fwd_cap : std::numeric_limits<capacity_type>::max();

    //    1-1: The source tree.
    for (p_node = p_s_start_node; ; ) {

        if ((p_fwd_arc = p_node->p_parent_arc) == TERMINAL) break;

        // Forward arc.
        if (!(p_node->tag & PARENT_REV)) {
            if (bottle_neck_cap > p_fwd_arc->rev_cap)
                bottle_neck_cap = p_fwd_arc->rev_cap;
            p_node = neighbor_node_fwd(p_node, p_fwd_arc->shift);
        }
        else {
            p_node = neighbor_node_rev(p_node, p_fwd_arc->shift);
        }

    } // for 

    if (bottle_neck_cap > p_node->cap)
        bottle_neck_cap = p_node->cap;
    if (bottle_neck_cap == 0) return;

    //    1-2: The sink tree.
    for (p_node = p_t_start_node; ; ) {

        if ((p_fwd_arc = p_node->p_parent_arc) == TERMINAL) break;
        
        // Reverse arc.
        if (p_node->tag & PARENT_REV) {
            if (bottle_neck_cap > p_fwd_arc->rev_cap)
                bottle_neck_cap = p_fwd_arc->rev_cap;
            p_node = neighbor_node_rev(p_node, p_fwd_arc->shift);
        }
        else {
            p_node = neighbor_node_fwd(p_node, p_fwd_arc->shift);               
        }

    } // for

    if (bottle_neck_cap > -p_node->cap)
        bottle_neck_cap = -p_node->cap;
    if (bottle_neck_cap == 0) return;


    // STEP 2: Augment
    if (0 != p_mid_rev_cap) *p_mid_rev_cap += bottle_neck_cap;
    if (0 != p_mid_fwd_cap) *p_mid_fwd_cap -= bottle_neck_cap;

    //    2-1: The source tree
    for (p_node = p_s_start_node; ; ) {

        if ((p_fwd_arc = p_node->p_parent_arc) == TERMINAL) break;

        if (!(p_node->tag & PARENT_REV)) {

            p_fwd_arc->rev_cap -= bottle_neck_cap;

            if (!p_fwd_arc->rev_cap) {
                // Add p_node to the orphan queue.
                p_node->p_parent_arc = (forward_arc_pointer)ORPHAN;
                m_orphan_nodes.push_back(p_node);
            }

            p_node = neighbor_node_fwd(p_node, p_fwd_arc->shift);
        }
        else {

            p_fwd_arc->rev_cap += bottle_neck_cap;

            p_node = neighbor_node_rev(p_node, p_fwd_arc->shift);
        }
    }

    p_node->cap -= bottle_neck_cap;
    if (0 == p_node->cap) {
        // Add p_node to the orphan queue.
        p_node->p_parent_arc = (forward_arc_pointer)ORPHAN;
        m_orphan_nodes.push_back(p_node);
    }

    //    2-2: The sink tree.
    for (p_node = p_t_start_node; ; ) {

        if ((p_fwd_arc = p_node->p_parent_arc) == TERMINAL) break;

        if (p_node->tag & PARENT_REV) {

            p_fwd_arc->rev_cap -= bottle_neck_cap;

            if (!p_fwd_arc->rev_cap) {
                p_node->p_parent_arc = (forward_arc_pointer)ORPHAN;
                m_orphan_nodes.push_back(p_node);
            }

            p_node = neighbor_node_rev(p_node, p_fwd_arc->shift);
        }
        else {
            
            p_fwd_arc->rev_cap += bottle_neck_cap;

            p_node = neighbor_node_fwd(p_node, p_fwd_arc->shift);
        }
    }

    p_node->cap += bottle_neck_cap;
    if (!p_node->cap) {

        // Add p_node to the orphan queue.
        p_node->p_parent_arc = (forward_arc_pointer)ORPHAN;
        m_orphan_nodes.push_back(p_node);
    }

    m_flow += bottle_neck_cap;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_fs_maxflow<_Cap, _Tg>::maxflow_adopt_source_orphan(
                                node_pointer p_orphan
                                )
{
    static const size_type
        DIST_MAX = std::numeric_limits<size_type>::max();

    node_pointer        p_node;
    reverse_arc_pointer p_first_in_arc, p_last_in_arc, p_rev_arc;
    forward_arc_pointer p_first_out_arc, p_last_out_arc, p_fwd_arc;
    forward_arc_pointer p_fwd_arc1, p_min_fwd_arc = 0;
    size_type           dist, min_dist = DIST_MAX;
    bool                parent_is_rev = false;


    if (p_orphan->tag & _Base::SPECIAL_OUT) {
        p_first_out_arc = p_orphan->p_first_out_arc + 1;
        p_last_out_arc
            = (forward_arc_pointer)(p_orphan->p_first_out_arc->shift);
    }
    else {
        p_first_out_arc = p_orphan->p_first_out_arc;
        p_last_out_arc  = (p_orphan + 1)->p_first_out_arc;
    }

    if (p_orphan->tag & _Base::SPECIAL_IN) {
        p_first_in_arc  = p_orphan->p_first_in_arc + 1;
        p_last_in_arc
            = (reverse_arc_pointer)(p_orphan->p_first_in_arc->p_fwd);
    }
    else {
        p_first_in_arc  = p_orphan->p_first_in_arc;
        p_last_in_arc   = (p_orphan + 1)->p_first_in_arc;
    }


    // Try to find new parent for the orphan.
    for (p_fwd_arc = p_first_out_arc;
         p_fwd_arc < p_last_out_arc;
         ++p_fwd_arc) {

        if (0 != p_fwd_arc->rev_cap) { // if 1

            p_node = neighbor_node_fwd(p_orphan, p_fwd_arc->shift);

            if (!(p_node->tag & IS_SINK)
            && (0 != p_node->p_parent_arc)) { // if 2

                dist = 0;

                while (true) {

                    if (p_node->dist_id == m_dist_id) {
                        dist += p_node->dist;
                        break;
                    }
                    ++dist;

                    p_fwd_arc1 = p_node->p_parent_arc;

                    if (p_fwd_arc1 == TERMINAL) {
                        p_node->dist_id = m_dist_id;
                        p_node->dist    = 1;
                        break;
                    } 
                    
                    if (p_fwd_arc1 == ORPHAN) {
                        dist = DIST_MAX;
                        break;
                    }

                    // Trace back to the source.
                    p_node = (!(p_node->tag & PARENT_REV)) ? 
                        neighbor_node_fwd(p_node, p_fwd_arc1->shift) : 
                            neighbor_node_rev(p_node, p_fwd_arc1->shift);

                } // while (true)

                if (dist < DIST_MAX) {

                    // Save minimum distance node so far.
                    if (dist < min_dist) {
                        parent_is_rev = false;
                        p_min_fwd_arc = p_fwd_arc;
                        min_dist      = dist;
                    }

                    // Set distance along the path.
                    for (p_node
                            = neighbor_node_fwd(p_orphan, p_fwd_arc->shift); 
                         p_node->dist_id != m_dist_id; ) {

                        p_node->dist_id = m_dist_id;
                        p_node->dist    = dist--;

                        p_fwd_arc1       = p_node->p_parent_arc;
                        
                        // Trace back to the source.
                        p_node = (!(p_node->tag & PARENT_REV)) ? 
                            neighbor_node_fwd(p_node, p_fwd_arc1->shift) : 
                                neighbor_node_rev(p_node, p_fwd_arc1->shift);

                    }
                } // if (dist < DIST_MAX)
            } // if 2
        } // if 1
    } // for

    for (p_rev_arc = p_first_in_arc;
         p_rev_arc < p_last_in_arc;
         ++p_rev_arc) {

        p_fwd_arc = p_rev_arc->p_fwd;

        p_node = neighbor_node_rev(p_orphan, p_fwd_arc->shift);

        if (!(p_node->tag & IS_SINK)
        && (0 != p_node->p_parent_arc)) { // if 4

            dist = 0;

            while (true) {

                if (p_node->dist_id == m_dist_id) {
                    dist += p_node->dist;
                    break;
                }
                ++dist;

                p_fwd_arc1 = p_node->p_parent_arc;

                if (p_fwd_arc1 == TERMINAL) {
                    p_node->dist_id = m_dist_id;
                    p_node->dist = 1;
                    break;
                } 
                
                if (p_fwd_arc1 == ORPHAN) {
                    dist = DIST_MAX;
                    break;
                }

                // Trace back to the source.
                p_node = (!(p_node->tag & PARENT_REV)) ? 
                    neighbor_node_fwd(p_node, p_fwd_arc1->shift) : 
                        neighbor_node_rev(p_node, p_fwd_arc1->shift);

            } // while (true)

            if (dist < DIST_MAX) {

                // Save minimum distance node so far.
                if (dist < min_dist) {
                    parent_is_rev = true;
                    p_min_fwd_arc = p_fwd_arc;
                    min_dist      = dist;
                }

                // Set distance along the path.
                for (p_node
                        = neighbor_node_rev(p_orphan, p_fwd_arc->shift); 
                     p_node->dist_id != m_dist_id; ) {

                    p_node->dist_id = m_dist_id;
                    p_node->dist    = dist--;

                    p_fwd_arc1       = p_node->p_parent_arc;
                    
                    // Trace back to the source.
                    p_node = (!(p_node->tag & PARENT_REV)) ? 
                        neighbor_node_fwd(p_node, p_fwd_arc1->shift) : 
                            neighbor_node_rev(p_node, p_fwd_arc1->shift);
                }
            } // if (dist < DIST_MAX)
        } // if 4

    } // for

    if (0 != (p_orphan->p_parent_arc = p_min_fwd_arc)) { // Found parent.

        if (!parent_is_rev) p_orphan->tag &= ~PARENT_REV;
        else                p_orphan->tag |=  PARENT_REV;

        p_orphan->dist_id  = m_dist_id;
        p_orphan->dist     = min_dist + 1;
    }
    else { // Parent not found.

        p_orphan->dist_id = 0;

        // Process neighbors.
        for (p_fwd_arc = p_first_out_arc;
             p_fwd_arc < p_last_out_arc;
             ++p_fwd_arc) {

            p_node = neighbor_node_fwd(p_orphan, p_fwd_arc->shift);

            if (!(p_node->tag & IS_SINK)
            && (0 != (p_fwd_arc1 = p_node->p_parent_arc))) {

                if (0 != p_fwd_arc->rev_cap)
                    activate(p_node);

                if ((p_fwd_arc1 != ORPHAN)
                && (p_fwd_arc1 != TERMINAL)
                && (p_node->tag & PARENT_REV)
                && (neighbor_node_rev(p_node, p_fwd_arc1->shift)
                    == p_orphan)) {

                    p_node->p_parent_arc = (forward_arc_pointer)ORPHAN;
                    m_orphan_nodes.push_back(p_node);
                }
            } // if
        } // for

        for (p_rev_arc = p_first_in_arc;
             p_rev_arc < p_last_in_arc;
             ++p_rev_arc) {

            p_fwd_arc = p_rev_arc->p_fwd;
            p_node   = neighbor_node_rev(p_orphan, p_fwd_arc->shift);
            
            if (!(p_node->tag & IS_SINK)
            && (0 != (p_fwd_arc1 = p_node->p_parent_arc))) {

                activate(p_node);

                if ((p_fwd_arc1 != ORPHAN)
                && (p_fwd_arc1 != TERMINAL)
                && (!(p_node->tag & PARENT_REV))
                && (neighbor_node_fwd(p_node, p_fwd_arc1->shift)
                    == p_orphan)) {

                    p_node->p_parent_arc = (forward_arc_pointer)ORPHAN;
                    m_orphan_nodes.push_back(p_node);
                }
            } // if
        } // for
    } // else
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
optnet_fs_maxflow<_Cap, _Tg>::maxflow_adopt_sink_orphan(
                                node_pointer p_orphan
                                )
{
    static const size_type
        DIST_MAX = std::numeric_limits<size_type>::max();

    node_pointer        p_node;
    reverse_arc_pointer p_first_in_arc, p_last_in_arc, p_rev_arc;
    forward_arc_pointer p_first_out_arc, p_last_out_arc, p_fwd_arc;
    forward_arc_pointer p_fwd_arc1, p_min_fwd_arc = 0;
    size_t              dist, min_dist = DIST_MAX;
    bool                parent_is_rev = false;

    if (p_orphan->tag & _Base::SPECIAL_OUT) {
        p_first_out_arc = p_orphan->p_first_out_arc + 1;
        p_last_out_arc 
            = (forward_arc_pointer)(p_orphan->p_first_out_arc->shift);
    }
    else {
        p_first_out_arc = p_orphan->p_first_out_arc;
        p_last_out_arc  = (p_orphan + 1)->p_first_out_arc;
    }

    if (p_orphan->tag & _Base::SPECIAL_IN) {
        p_first_in_arc  = p_orphan->p_first_in_arc + 1;
        p_last_in_arc 
            = (reverse_arc_pointer)(p_orphan->p_first_in_arc->p_fwd);
    }
    else {
        p_first_in_arc  = p_orphan->p_first_in_arc;
        p_last_in_arc   = (p_orphan + 1)->p_first_in_arc;
    }

    // Try to find new parent for the orphan.
    for (p_fwd_arc = p_first_out_arc;
         p_fwd_arc < p_last_out_arc;
         ++p_fwd_arc) {

        p_node = neighbor_node_fwd(p_orphan, p_fwd_arc->shift);

        if ((p_node->tag & IS_SINK) && (p_node->p_parent_arc)) { // if 2

            dist = 0;

            while (true) {

                if (p_node->dist_id == m_dist_id) {
                    dist += p_node->dist;
                    break;
                }
                ++dist;

                p_fwd_arc1 = p_node->p_parent_arc;

                if (p_fwd_arc1 == TERMINAL) {
                    p_node->dist_id = m_dist_id;
                    p_node->dist = 1;
                    break;
                } 
                
                if (p_fwd_arc1 == ORPHAN) {
                    dist = DIST_MAX;
                    break;
                }

                // Trace back to the source.
                p_node = (!(p_node->tag & PARENT_REV)) ? 
                    neighbor_node_fwd(p_node, p_fwd_arc1->shift) : 
                        neighbor_node_rev(p_node, p_fwd_arc1->shift);

            } // while (true)

            if (dist < DIST_MAX) {

                // Save minimum distance node so far.
                if (dist < min_dist) {
                    parent_is_rev = false;
                    p_min_fwd_arc = p_fwd_arc;
                    min_dist      = dist;
                }

                // Set distance along the path.
                for (p_node
                        = neighbor_node_fwd(p_orphan, p_fwd_arc->shift); 
                     p_node->dist_id != m_dist_id; ) {

                    p_node->dist_id = m_dist_id;
                    p_node->dist    = dist--;

                    p_fwd_arc1       = p_node->p_parent_arc;
                    
                    // Trace back to the source.
                    p_node = (!(p_node->tag & PARENT_REV)) ? 
                        neighbor_node_fwd(p_node, p_fwd_arc1->shift) : 
                            neighbor_node_rev(p_node, p_fwd_arc1->shift);

                }
            } // if (dist < DIST_MAX)
        } // if 2
    } // for

    for (p_rev_arc = p_first_in_arc;
         p_rev_arc < p_last_in_arc;
         ++p_rev_arc) {

        p_fwd_arc = p_rev_arc->p_fwd;

        if (0 != p_fwd_arc->rev_cap) { // if 3
            
            p_node = neighbor_node_rev(p_orphan, p_fwd_arc->shift);

            if ((p_node->tag & IS_SINK) && (p_node->p_parent_arc)) { // if 4

                dist = 0;

                while (true) {

                    if (p_node->dist_id == m_dist_id) {
                        dist += p_node->dist;
                        break;
                    }
                    ++dist;

                    p_fwd_arc1 = p_node->p_parent_arc;

                    if (p_fwd_arc1 == TERMINAL) {
                        p_node->dist_id = m_dist_id;
                        p_node->dist = 1;
                        break;
                    } 
                
                    if (p_fwd_arc1 == ORPHAN) {
                        dist = DIST_MAX;
                        break;
                    }

                    // Trace back to the source.
                    p_node = (!(p_node->tag & PARENT_REV)) ? 
                        neighbor_node_fwd(p_node, p_fwd_arc1->shift) : 
                            neighbor_node_rev(p_node, p_fwd_arc1->shift);

                } // while (true)

                if (dist < DIST_MAX) {

                    // Save minimum distance node so far.
                    if (dist < min_dist) {
                        parent_is_rev = true;
                        p_min_fwd_arc = p_fwd_arc;
                        min_dist      = dist;
                    }

                    // Set distance along the path.
                    for (p_node = neighbor_node_rev(p_orphan, p_fwd_arc->shift); 
                            p_node->dist_id != m_dist_id; ) {

                        p_node->dist_id = m_dist_id;
                        p_node->dist    = dist--;

                        p_fwd_arc1       = p_node->p_parent_arc;
                    
                        // Trace back to the source.
                        p_node = (!(p_node->tag & PARENT_REV)) ? 
                            neighbor_node_fwd(p_node, p_fwd_arc1->shift) : 
                                neighbor_node_rev(p_node, p_fwd_arc1->shift);
                    }
                } // if (dist < DIST_MAX)
            } // if 4
        } // if 3
    } // for

    if (0 != (p_orphan->p_parent_arc = p_min_fwd_arc)) { // Found parent.

        if (!parent_is_rev) p_orphan->tag &= ~PARENT_REV;
        else                p_orphan->tag |=  PARENT_REV;

        p_orphan->dist_id = m_dist_id;
        p_orphan->dist    = min_dist + 1;
    }
    else { // Parent not found.

        p_orphan->dist_id = 0;

        // Process neighbors.
        for (p_fwd_arc = p_first_out_arc;
             p_fwd_arc < p_last_out_arc;
             ++p_fwd_arc) {

            p_node = neighbor_node_fwd(p_orphan, p_fwd_arc->shift);

            if ((p_node->tag & IS_SINK)
            && (0 != (p_fwd_arc1 = p_node->p_parent_arc))) {

                activate(p_node);

                if ((p_fwd_arc1 != ORPHAN)
                && (p_fwd_arc1 != TERMINAL)
                && (p_node->tag & PARENT_REV)
                && (neighbor_node_rev(p_node, p_fwd_arc1->shift)
                        == p_orphan)) {

                    p_node->p_parent_arc = (forward_arc_pointer)ORPHAN;
                    m_orphan_nodes.push_back(p_node);
                }
            } // if
        } // for

        for (p_rev_arc = p_first_in_arc;
             p_rev_arc < p_last_in_arc;
             ++p_rev_arc) {

            p_fwd_arc = p_rev_arc->p_fwd;
            p_node   = neighbor_node_rev(p_orphan, p_fwd_arc->shift);
            
            if ((p_node->tag & IS_SINK)
            && (0 != (p_fwd_arc1 = p_node->p_parent_arc))) {

                if (0 != p_fwd_arc->rev_cap)
                    activate(p_node);

                if ((p_fwd_arc1 != ORPHAN)
                && (p_fwd_arc1 != TERMINAL)
                && (!(p_node->tag & PARENT_REV))
                && (neighbor_node_fwd(p_node, p_fwd_arc1->shift)
                        == p_orphan)) {

                    p_node->p_parent_arc = (forward_arc_pointer)ORPHAN;
                    m_orphan_nodes.push_back(p_node);
                }
            } // if
        } // for
    } // else
}

//
// Constants
//
template<typename _Cap, typename _Tg>
    const unsigned char optnet_fs_maxflow<_Cap, _Tg>::PARENT_REV = 0x10;

template<typename _Cap, typename _Tg>
    const unsigned char optnet_fs_maxflow<_Cap, _Tg>::IS_ACTIVE  = 0x04;

template<typename _Cap, typename _Tg>
    const unsigned char optnet_fs_maxflow<_Cap, _Tg>::IS_SINK    = 0x01;

template<typename _Cap, typename _Tg>
    typename optnet_fs_maxflow<_Cap, _Tg>::forward_arc_const_pointer
        optnet_fs_maxflow<_Cap, _Tg>::TERMINAL = 
            (typename optnet_fs_maxflow<_Cap, _Tg>::forward_arc_const_pointer)(1);

template<typename _Cap, typename _Tg>
    typename optnet_fs_maxflow<_Cap, _Tg>::forward_arc_const_pointer
        optnet_fs_maxflow<_Cap, _Tg>::ORPHAN = 
            (typename optnet_fs_maxflow<_Cap, _Tg>::forward_arc_const_pointer)(2);

} // namespace

#endif
