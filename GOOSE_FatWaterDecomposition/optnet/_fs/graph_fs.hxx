/*
 ==========================================================================
 |   
 |   $Id: graph_fs.hxx 2137 2007-07-02 03:26:31Z kangli $
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
     
     This file implements a graph data structure using the forward-star
     representation.

 ==========================================================================
 */

#ifndef ___GRAPH_FS_HXX___
#   define ___GRAPH_FS_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/config.h>

#   ifdef __OPTNET_SUPPORT_ROI__
#       include <optnet/define.h> // for roi_node definition
#   endif

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/chunk_list.hxx>
#   include <optnet/_base/graph_traits.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class graph_fs
///  @brief A graph class using the forward-star representation scheme.
///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg = net_f_xy>
class graph_fs
{
    struct fwd_arc;
    struct rev_arc;

    struct node
    {
        fwd_arc*        p_first_out_arc; // Pointer to the first out arc.
        rev_arc*        p_first_in_arc;  // Pointer to the first in arc.
        fwd_arc*        p_parent_arc;    // Pointer to the parent arc.
        size_t          dist_id;
        size_t          dist;
        _Cap            cap;
        unsigned char   tag;
    };

    // Forward arc struct.
    struct fwd_arc
    {
        ptrdiff_t       shift;          // Tail node + shift = head node.
        _Cap            rev_cap;        // Reverse residual capacity.
    };

    // Reverse arc struct.
    struct rev_arc
    {
        fwd_arc*        p_fwd;          // Pointer to the forward arc.
    };


public:

    typedef array<node, _Tg>                    node_container;
    typedef chunk_list<fwd_arc>                 forward_arc_container;
    typedef chunk_list<rev_arc>                 reverse_arc_container;

    typedef graph_traits<
        _Cap, 
        node_container, 
        forward_arc_container, 
        reverse_arc_container>                  _Traits;

    typedef typename
        _Traits::node_type                      node_type;
    typedef typename 
        _Traits::node_reference                 node_reference;
    typedef typename 
        _Traits::node_const_reference           node_const_reference;
    typedef typename 
        _Traits::node_iterator                  node_iterator;
    typedef typename 
        _Traits::node_const_iterator            node_const_iterator;
    typedef typename 
        _Traits::node_pointer                   node_pointer;
    typedef typename 
        _Traits::node_const_pointer             node_const_pointer;

    typedef typename 
        _Traits::forward_arc_type               forward_arc_type;
    typedef typename 
        _Traits::forward_arc_reference          forward_arc_reference;
    typedef typename 
        _Traits::forward_arc_const_reference    forward_arc_const_reference;
    typedef typename 
        _Traits::forward_arc_iterator           forward_arc_iterator;
    typedef typename 
        _Traits::forward_arc_const_iterator     forward_arc_const_iterator;
    typedef typename 
        _Traits::forward_arc_pointer            forward_arc_pointer;
    typedef typename 
        _Traits::forward_arc_const_pointer      forward_arc_const_pointer;

    typedef typename 
        _Traits::reverse_arc_type               reverse_arc_type;
    typedef typename 
        _Traits::reverse_arc_reference          reverse_arc_reference;
    typedef typename 
        _Traits::reverse_arc_const_reference    reverse_arc_const_reference;
    typedef typename 
        _Traits::reverse_arc_iterator           reverse_arc_iterator;
    typedef typename 
        _Traits::reverse_arc_const_iterator     reverse_arc_const_iterator;
    typedef typename 
        _Traits::reverse_arc_pointer            reverse_arc_pointer;
    typedef typename 
        _Traits::reverse_arc_const_pointer      reverse_arc_const_pointer;

    typedef typename _Traits::capacity_type     capacity_type;
    typedef typename _Traits::difference_type   difference_type;
    typedef typename _Traits::size_type         size_type;

    // ROI definitions.
#   ifdef __OPTNET_SUPPORT_ROI__
    typedef roi_node_t                          roi_node_type;
    typedef array_base<roi_node_t>              roi_base_type;
    typedef array_ref<roi_node_t>               roi_ref_type;
    typedef array<roi_node_t>                   roi_type;
#   endif

    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    graph_fs();

    ///////////////////////////////////////////////////////////////////////
    ///  Construct a graph_fs object of the given size.
    ///
    ///  @param  s0   The size of the first  dimension of the graph. 
    ///  @param  s1   The size of the second dimension of the graph. 
    ///  @param  s2   The size of the third  dimension of the graph. 
    ///  @param  s3   The size of the fourth dimension of the graph
    ///               (default: 1).
    ///  @param  s4   The size of the fifth  dimension of the graph
    ///               (default: 1).
    ///
    ///////////////////////////////////////////////////////////////////////
    graph_fs(size_type s0, 
             size_type s1, 
             size_type s2, 
             size_type s3 = 1, 
             size_type s4 = 1
             );

    ///////////////////////////////////////////////////////////////////////
    /// Default destructor.
    ///////////////////////////////////////////////////////////////////////
    virtual ~graph_fs() {}

    ///////////////////////////////////////////////////////////////////////
    ///  Create a graph of the given size.
    ///
    ///  @param  s0   The size of the first  dimension of the graph. 
    ///  @param  s1   The size of the second dimension of the graph. 
    ///  @param  s2   The size of the third  dimension of the graph. 
    ///  @param  s3   The size of the fourth dimension of the graph
    ///               (default: 1).
    ///  @param  s4   The size of the fifth  dimension of the graph
    ///               (default: 1).
    ///
    ///  @returns Returns true if the graph is successfully created,
    ///           false otherwise.
    ///
    ///  @remarks The created graph can have up to five dimensions.
    ///
    ///////////////////////////////////////////////////////////////////////
    virtual bool create(size_type s0, 
                        size_type s1, 
                        size_type s2, 
                        size_type s3 = 1,
                        size_type s4 = 1
                        );

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
    ///  Returns the size of the fourth dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_3() const { return m_nodes.size_3(); }
    
    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of the fifth  dimension of the graph.
    ///////////////////////////////////////////////////////////////////////
    inline size_type size_4() const { return m_nodes.size_4(); }


#   ifdef __OPTNET_SUPPORT_ROI__

    ///////////////////////////////////////////////////////////////////////
    ///  Set the region of interest.
    ///
    ///  @param roi The region-of-interest mask array.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline bool set_roi(const roi_base_type& roi)
    {
        //FIXME: Only considered up to 4-D.
        if (roi.size_0() != size_0() || 
            roi.size_1() != size_1() ||
            roi.size_2() != size_3()
            )
            return false;
        
        m_proi = &roi;
        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns whether an ROI has been specified.
    ///////////////////////////////////////////////////////////////////////
    inline bool has_roi() const { return (m_proi != 0); }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns whether the given position is in the ROI.
    ///////////////////////////////////////////////////////////////////////
    inline bool in_roi(size_type i0,
                       size_type i1,
                       size_type i2,
                       size_type i3 = 0,
                       size_type i4 = 0 // <-- not used
                       ) const
    {
        OPTNET_UNUSED(i4);
        //FIXME: Only considered up to 4-D.
        return (
                !has_roi() || (
                    i2 >= (*m_proi)(i0, i1, i3).lower &&
                    i2 <  (*m_proi)(i0, i1, i3).upper
                )
            );
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns whether the speficied position is the lowest in the
    ///  column in the ROI.
    ///////////////////////////////////////////////////////////////////////
    inline bool is_lowest(size_type i0,
                          size_type i1,
                          size_type i2,
                          size_type i3 = 0,
                          size_type i4 = 0 // <-- not used
                          ) const
    {
        OPTNET_UNUSED(i4);
        assert(0 != m_proi);
        return ((*m_proi)(i0, i1, i3).lower == i2);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to an ROI node.
    ///////////////////////////////////////////////////////////////////////
    inline const roi_node_type& roi_node(size_type i0,
                                         size_type i1,
                                         size_type i3 = 0,
                                         size_type i4 = 0 // <-- not used
                                         ) const
    {
        OPTNET_UNUSED(i4);
        assert(0 != m_proi);
        return (*m_proi)(i0, i1, i3);
    }

#   endif // __OPTNET_SUPPORT_ROI__


    ///////////////////////////////////////////////////////////////////////
    ///  Add an arc connecting from the specified tail node to the
    ///  specified head node.
    ///
    ///  @param  tail_i0 The first  index of the tail node.
    ///  @param  tail_i1 The second index of the tail node.
    ///  @param  tail_i2 The third  index of the tail node.
    ///  @param  head_i0 The first  index of the head node.
    ///  @param  head_i1 The second index of the head node.
    ///  @param  head_i2 The third  index of the head node.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc(size_type tail_i0,
                        size_type tail_i1,
                        size_type tail_i2,
                        size_type head_i0,
                        size_type head_i1,
                        size_type head_i2
                        )
    {
#   ifdef __OPTNET_SUPPORT_ROI__
        if (!in_roi(tail_i0, tail_i1, tail_i2) ||
            !in_roi(head_i0, head_i1, head_i2)
            ) {
            //...
            return;
        }
#   endif

        node_pointer p_tail_node = &(m_nodes(tail_i0,
                                             tail_i1,
                                             tail_i2
                                             ));
        node_pointer p_head_node = &(m_nodes(head_i0,
                                             head_i1,
                                             head_i2
                                             ));

        add_arc_helper(p_tail_node, p_head_node);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add an arc connecting from the specified tail node to the
    ///  specified head node.
    ///
    ///  @param  tail_i0 The first  index of the tail node.
    ///  @param  tail_i1 The second index of the tail node.
    ///  @param  tail_i2 The third  index of the tail node.
    ///  @param  tail_i3 The fourth index of the tail node.
    ///  @param  head_i0 The first  index of the head node.
    ///  @param  head_i1 The second index of the head node.
    ///  @param  head_i2 The third  index of the head node.
    ///  @param  head_i3 The fourth index of the head node.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc(size_type tail_i0,
                        size_type tail_i1,
                        size_type tail_i2,
                        size_type tail_i3,
                        size_type head_i0,
                        size_type head_i1,
                        size_type head_i2,
                        size_type head_i3
                        )
    {
#   ifdef __OPTNET_SUPPORT_ROI__
        if (!in_roi(tail_i0, tail_i1, tail_i2, tail_i3) ||
            !in_roi(head_i0, head_i1, head_i2, head_i3)
            ) {
            //...
            return;
        }
#   endif

        node_pointer p_tail_node = &(m_nodes(tail_i0,
                                             tail_i1,
                                             tail_i2,
                                             tail_i3
                                             ));
        node_pointer p_head_node = &(m_nodes(head_i0,
                                             head_i1,
                                             head_i2,
                                             head_i3
                                             ));

        add_arc_helper(p_tail_node, p_head_node);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add an arc connecting from the specified tail node to the
    ///  specified head node.
    ///
    ///  @param  tail_i0 The first  index of the tail node.
    ///  @param  tail_i1 The second index of the tail node.
    ///  @param  tail_i2 The third  index of the tail node.
    ///  @param  tail_i3 The fourth index of the tail node.
    ///  @param  tail_i4 The fifth  index of the tail node.
    ///  @param  head_i0 The first  index of the head node.
    ///  @param  head_i1 The second index of the head node.
    ///  @param  head_i2 The third  index of the head node.
    ///  @param  head_i3 The fourth index of the head node.
    ///  @param  head_i4 The fifth  index of the head node.
    ///
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc(size_type tail_i0,
                        size_type tail_i1,
                        size_type tail_i2,
                        size_type tail_i3,
                        size_type tail_i4,
                        size_type head_i0,
                        size_type head_i1,
                        size_type head_i2,
                        size_type head_i3,
                        size_type head_i4
                        )
    {
#   ifdef __OPTNET_SUPPORT_ROI__
        if (!in_roi(tail_i0, tail_i1, tail_i2, tail_i3, tail_i4) ||
            !in_roi(head_i0, head_i1, head_i2, head_i3, head_i4)
            ) {
            //...
            return;
        }
#   endif

        node_pointer p_tail_node = &(m_nodes(tail_i0,
                                             tail_i1,
                                             tail_i2,
                                             tail_i3,
                                             tail_i4
                                             ));
        node_pointer p_head_node = &(m_nodes(head_i0,
                                             head_i1,
                                             head_i2,
                                             head_i3,
                                             head_i4
                                             ));

        add_arc_helper(p_tail_node, p_head_node);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Remove all arcs that were added.
    ///////////////////////////////////////////////////////////////////////
    inline void clear_arcs()
    {
        if (m_fwd_arcs.count() > 0 || m_rev_arcs.count() > 0) {

            m_fwd_arcs.clear();
            m_rev_arcs.clear();

            for (node_iterator it = m_nodes.begin();
                 it != m_nodes.end();
                 ++it) {

                it->p_first_out_arc
                    = reinterpret_cast<fwd_arc*>(0);
                it->p_first_in_arc
                    = reinterpret_cast<rev_arc*>(0);
            }
        }
        
        // Set "prepared" flag to false.
        // Call prepare() to make the graph usable.
        m_prepared = false;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Converts the arcs added by 'add_arc()' calls to a forward-star
    ///  representaion, i.e., the arcs are sorted in the order of their
    ///  tail nodes.
    ///
    ///  @remarks After calling this function, the graph is ready for
    ///           use, but at the mean time, one cannot call the
    ///           function 'add_arc()' any more.
    ///////////////////////////////////////////////////////////////////////
    void prepare();


protected:

    ///////////////////////////////////////////////////////////////////////
    ///  This is an helper function to the 'add_arc()' calls.
    ///////////////////////////////////////////////////////////////////////
    inline void add_arc_helper(node_pointer p_tail_node,
                               node_pointer p_head_node
                               )
    {
        assert(m_prepared == false);

        // Allocate new forward and reverse arcs and append them to the
        // corresponding arc lists.
        rev_arc* p_rev_arc   = m_rev_arcs.append();
        fwd_arc* p_fwd_arc   = m_fwd_arcs.append();

        // Temporarily store the tail and head nodes in the p_fwd and
        // shift fields of the reverse and forward arcs.
        p_rev_arc->p_fwd     = reinterpret_cast<fwd_arc *>(p_tail_node);
        p_fwd_arc->shift     = reinterpret_cast<ptrdiff_t>(p_head_node);

        // The reverse residual capacity is initialized to be zero.
        p_fwd_arc->rev_cap   = static_cast<capacity_type>(0);
        
        // Temporarily store the number of outgoing and incoming arcs
        // in the p_first_out_arc and p_first_in_arc fields.
        p_tail_node->p_first_out_arc = reinterpret_cast<fwd_arc*>(
            reinterpret_cast<ptrdiff_t>
                (p_tail_node->p_first_out_arc) + 1
            );
        p_head_node->p_first_in_arc = reinterpret_cast<rev_arc*>(
            reinterpret_cast<ptrdiff_t>
                (p_head_node->p_first_in_arc) + 1
            );
    }

    // Constants (Initialized in graph_fs.cxx)
    static const unsigned char  SPECIAL_OUT;
    static const unsigned char  SPECIAL_IN;

#   ifdef __OPTNET_SUPPORT_ROI__
    const roi_base_type*        m_proi;
#   endif // __OPTNET_SUPPORT_ROI__

    bool                        m_prepared;
    node_container              m_nodes;
    forward_arc_container       m_fwd_arcs;
    reverse_arc_container       m_rev_arcs;
};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_fs/graph_fs.cxx>
#   endif

#endif
