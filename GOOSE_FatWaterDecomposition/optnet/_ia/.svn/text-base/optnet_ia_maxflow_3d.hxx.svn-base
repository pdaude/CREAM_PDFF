/*
 ==========================================================================
 |   
 |   $Id: optnet_ia_maxflow_3d.hxx 2137 2007-07-02 03:26:31Z kangli $
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

/*
 ==========================================================================
  - Purpose:

      This file implements Boykov--Kolmogorov's max-flow/min-cut algorithm
      on a graph whose arcs are implicitly represented.
      
      The implementation makes one or more simplifications, which are only
      applicable to the Optimal-Net algorithm, as listed below:

         * The residual capacity of non-st arcs is assumed to be +infinity
           and will not be decremented.


  - Reference(s):

    [1] Yuri Boykov and Vladimir Kolmogorov
        An Experimental Comparison of Min-Cut/Max-Flow Algorithms for 
            Energy Minimization in Vision
        IEEE Trans. on Pattern Analysis and Machine Intelligence, 2004
        URL: http://www.csd.uwo.ca/faculty/yuri/Abstracts/pami04-abs.html
 ==========================================================================
 */

#ifndef ___OPTNET_IA_MAXFLOW_3D_HXX___
#   define ___OPTNET_IA_MAXFLOW_3D_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4127)
#   endif

#   include <optnet/define.h>
#   include <optnet/_base/graph_traits.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/array.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <algorithm>
#   include <queue>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_ia_maxflow_3d
///  @brief 3-D implicit-arc implementation of the Boykov-Kolmogorov
///         max-flow algorithm.
///  @see   optnet_ia_3d
///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg = net_f_xy>
class optnet_ia_maxflow_3d
{
    struct node
    {
        node*   p_parent;       // Parent node pointer.
        size_t  dist_id;        // Distance ID.
        size_t  dist;           // Distance value.
        long    tag;            // Tag for marking properties of node.
        _Cap    res_cap[5];     // Reversed s-t arc capacity.
        _Cap    cap;            // S-t arc capacity.
    };

    // These are used internally.
    static const node*                              TERMINAL;
    static const node*                              ORPHAN;

    static const long                               FWD_ARC_MASK;
    static const long                               REV_ARC_MASK;
    static const long                               PARENT_MASK;
    static const long                               PARENT_REV;
    static const long                               IS_ACTIVE;
    static const long                               IS_SINK;


public:


    typedef array<node, _Tg>                        node_container;

    typedef graph_traits<_Cap, node_container>      _Traits;

    typedef roi_node_t                              roi_node_type;
    typedef array_base<roi_node_t>                  roi_base_type;
    typedef array_ref<roi_node_t>                   roi_ref_type;
    typedef array<roi_node_t>                       roi_type;

    typedef typename _Traits::node_type             node_type;
    typedef typename _Traits::node_reference        node_reference;
    typedef typename _Traits::node_const_reference  node_const_reference;
    typedef typename _Traits::node_iterator         node_iterator;
    typedef typename _Traits::node_const_iterator   node_const_iterator;
    typedef typename _Traits::node_pointer          node_pointer;
    typedef typename _Traits::node_const_pointer    node_const_pointer;
    typedef typename _Traits::capacity_type         capacity_type;
    typedef typename _Traits::difference_type       difference_type;
    typedef typename _Traits::size_type             size_type;

    typedef std::deque<node_pointer>                node_queue;


    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    optnet_ia_maxflow_3d();

    ///////////////////////////////////////////////////////////////////////
    ///  Construct a optnet_ia_maxflow_3d object with the underlying graph
    ///  being created according to the given size information.
    ///
    ///  @param  s0   The size of the first  dimension of the graph. 
    ///  @param  s1   The size of the second dimension of the graph. 
    ///  @param  s2   The size of the third  dimension of the graph. 
    ///
    ///////////////////////////////////////////////////////////////////////
    optnet_ia_maxflow_3d(size_type s0, size_type s1, size_type s2);
    
    ///////////////////////////////////////////////////////////////////////
    ///  Construct a optnet_pr_maxflow object with the underlying graph
    ///  being created according to the given size information.
    ///
    ///  @param  s0   The size of the first  dimension of the graph. 
    ///  @param  s1   The size of the second dimension of the graph. 
    ///  @param  s2   The size of the third  dimension of the graph. 
    ///
    ///////////////////////////////////////////////////////////////////////
    bool create(size_type s0, size_type s1, size_type s2);

    ///////////////////////////////////////////////////////////////////////
    ///  Solve the maximum-flow/minimum s-t cut problem.
    ///
    ///  @returns The maximum flow value.
    ///////////////////////////////////////////////////////////////////////
    capacity_type solve();

    ///////////////////////////////////////////////////////////////////////
    ///  Set the smoothness constraints.
    ///
    ///  @param smooth0 The first smoothness parameter.
    ///  @param smooth1 The second smoothness parameter.
    ///  @param circle0 Enabling/disabling circle graph construction.
    ///  @param circle1 Enabling/disabling circle graph construction.
    ///
    ///  @remarks The actual meanings of the smoothness parameters depend
    ///           on the orientation setting of the net.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool set_params(int  smooth0 = 1,
                           int  smooth1 = 1,
                           bool circle0 = false,
                           bool circle1 = false
                           )
    {
        if (smooth0 < 0 || smooth1 < 0) return false;

        m_circle[0] = circle0;
        m_circle[1] = circle1;

        m_smooth[0] = smooth0;
        m_smooth[1] = smooth1;

        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Set the region of interest.
    ///
    ///  @param roi The region-of-interest mask array.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool set_roi(const roi_base_type& roi)
    {
        if (roi.size_0() != size_0() || 
            roi.size_1() != size_1())
            return false;
        
        m_proi = &roi;
        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add arc(s) connecting a node to the source and/or the sink node.
    ///
    ///  @param  s    The capacity of the arc from the source node.
    ///  @param  t    The capacity of the arc to the sink node.
    ///  @param  i0   The first  index of the node. 
    ///  @param  i1   The second index of the node. 
    ///  @param  i2   The third  index of the node. 
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_st_arc(capacity_type s,
                           capacity_type t,
                           size_type     i0,
                           size_type     i1,
                           size_type     i2
                           )
    {
        if (in_roi(i0, i1, i2)) {
            m_nodes(i0, i1, i2).cap = s - t;
            m_preflow += (s < t) ? s : t;
        }
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Determines if the given node is in the source set of the cut.
    ///
    ///  @param  i0   The first  index of the node. 
    ///  @param  i1   The second index of the node. 
    ///  @param  i2   The third  index of the node. 
    ///
    ///  @return Returns true if the given node is the source set,
    ///          false otherwise.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool in_source_set(size_type i0,
                              size_type i1,
                              size_type i2
                              ) const
    {
        return (0 == (m_nodes(i0, i1, i2).tag & IS_SINK)) &&
               (0 != m_nodes(i0, i1, i2).p_parent) &&
               in_roi(i0, i1, i2);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the total number of nodes in the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size()   const { return m_nodes.size();   }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the first  dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_0() const { return m_nodes.size_0(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the second dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_1() const { return m_nodes.size_1(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the third  dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_2() const { return m_nodes.size_2(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns whether an ROI has been specified.
    ///////////////////////////////////////////////////////////////////////
    inline bool has_roi() const { return (m_proi != 0); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns whether the given column has portions in the ROI.
    ///////////////////////////////////////////////////////////////////////
    inline bool in_roi(size_type i0, size_type i1) const
    {
        return (
                !has_roi() || (
                    0 != (*m_proi)(i0, i1, 0).lower &&
                    0 != (*m_proi)(i0, i1, 0).upper
                )
            );
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns whether the given position is in the ROI.
    ///////////////////////////////////////////////////////////////////////
    inline bool in_roi(size_type i0, size_type i1, size_type i2) const
    {
        return (
                !has_roi() || (
                    i2 >= (*m_proi)(i0, i1, 0).lower &&
                    i2 <  (*m_proi)(i0, i1, 0).upper
                )
            );
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns whether the speficied position is the lowest in the
    ///  column in the ROI.
    ///////////////////////////////////////////////////////////////////////
    inline bool is_lowest(size_type i0, size_type i1, size_type i2) const
    {
        assert(0 != m_proi);
        return ((*m_proi)(i0, i1, 0).lower == i2);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the index of the lowest node in column(i0, i1)
    ///////////////////////////////////////////////////////////////////////
    inline size_type lowest(size_type i0, size_type i1) const
    {
        if (0 != m_proi)
            return (*m_proi)(i0, i1, 0).lower;
        return 0;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Set the initial flow value.
    ///
    ///  @param flow The initial flow value.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void set_initial_flow(const capacity_type& flow)
    {
        m_preflow = flow;
    }


private:

    void init();

    void maxflow_init();
    void maxflow_augment(node_pointer   p_s_start_node,
                         node_pointer   p_t_start_node,
                         capacity_type* p_cap_mid,
                         capacity_type* p_mid_rev_cap);
    void maxflow_adopt_source_orphan(node_pointer p_orphan);
    void maxflow_adopt_sink_orphan(node_pointer p_orphan);

    inline void activate(node_pointer p_node)
    {
        if (0 == (p_node->tag & IS_ACTIVE)) {  // Not active yet.
            m_active_nodes.push_back(p_node);
            p_node->tag |= IS_ACTIVE;
        }
    }

    inline int fwd_arc_index(long tag)
    {
        return (int)((tag & FWD_ARC_MASK) >> 4);
    }

    inline int rev_arc_index(long tag)
    {
        return (int)((tag & REV_ARC_MASK) >> 8);
    }

    const roi_base_type*    m_proi;
    int                     m_max_arcs;
    int                     m_smooth[2];
    bool                    m_circle[2];
    difference_type         m_offset[10][5];
    long                    m_mask_in[5], m_mask_out[5];
    node_queue              m_active_nodes, m_orphan_nodes;
    capacity_type           m_preflow, m_flow;
    size_type               m_dist_id;
    node_container          m_nodes;
};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_ia/optnet_ia_maxflow_3d.cxx>
#   endif

#endif
