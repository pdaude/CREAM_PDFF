/*
 ==========================================================================
 |   
 |   $Id: graph_fs.cxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___GRAPH_FS_CXX___
#   define ___GRAPH_FS_CXX___

#   include <optnet/_fs/graph_fs.hxx>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
graph_fs<_Cap, _Tg>::graph_fs() :
#   ifdef __OPTNET_SUPPORT_ROI__
    m_proi(0),
#   endif // __OPTNET_SUPPORT_ROI__
    m_prepared(false)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
graph_fs<_Cap, _Tg>::graph_fs(size_type s0, 
                              size_type s1, 
                              size_type s2, 
                              size_type s3, 
                              size_type s4
                              ) :
#   ifdef __OPTNET_SUPPORT_ROI__    
    m_proi(0),
#   endif // __OPTNET_SUPPORT_ROI__
    m_prepared(false),
    m_nodes(s0, s1, s2, s3, s4)
{
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
bool 
graph_fs<_Cap, _Tg>::create(size_type s0, 
                            size_type s1, 
                            size_type s2, 
                            size_type s3,
                            size_type s4
                            )
{
    // Clears the ROI (Shouldn't do this!).
    // m_proi = 0;
    
    // Remove all arcs if any.
    clear_arcs();

    return m_nodes.create(s0, s1, s2, s3, s4);
}

///////////////////////////////////////////////////////////////////////////
//
// Convert arcs added by 'add_arc()' calls to a forward-star representaion,
// i.e., the arcs are sorted by the order of their tail nodes.
//
// Linear time algorithm.
// No or little additional memory will be allocated during the conversion.
// Post-condition: Arcs corresponding to the same node must be contiguous,
//                 i.e., be in one arc chunk.
//
///////////////////////////////////////////////////////////////////////////
template <typename _Cap, typename _Tg>
void
graph_fs<_Cap, _Tg>::prepare()
{
    typename forward_arc_container::chunk_pointer p_fwd_chunk;
    typename reverse_arc_container::chunk_pointer p_rev_chunk;
    typename forward_arc_container::list_type::iterator 
                it_fwd_chunk_ptr, it_fwd_chunk_ptr_end;
    typename reverse_arc_container::list_type::iterator 
                it_rev_chunk_ptr, it_rev_chunk_ptr_end,
                it_rev_chunk_ptr1;

    bool                rev_is_full, fwd_is_full;
    node_pointer        p_node, p_node1, p_tail_node, p_head_node;
    forward_arc_pointer p_fwd_arc, p_fwd_arc1;
    reverse_arc_pointer p_rev_arc, p_rev_arc1;

    if (m_prepared) return;

    ///////////////////////////////////////////////////////////////////////
    // Assertions:
    //   1) Arc container cannot be empty.
    //   2) # of forward arcs must match # of reverse arcs.
    assert(m_fwd_arcs.list().size() != 0);
    assert(m_fwd_arcs.list().size() == m_rev_arcs.list().size());

    it_fwd_chunk_ptr     = m_fwd_arcs.list().begin();
    it_rev_chunk_ptr     = m_rev_arcs.list().begin();
    it_rev_chunk_ptr1    = it_rev_chunk_ptr;

    it_fwd_chunk_ptr_end = m_fwd_arcs.list().end();
    it_rev_chunk_ptr_end = m_rev_arcs.list().end();

    p_fwd_chunk          = *it_fwd_chunk_ptr;
    p_rev_chunk          = *it_rev_chunk_ptr;

    p_fwd_arc            = p_fwd_chunk->begin();
    p_rev_arc            = p_rev_chunk->begin();
    p_rev_arc1           = p_rev_arc;

    // Flags marking the beginning of allocation of new storages.
    fwd_is_full          = false;
    rev_is_full          = false;

    ///////////////////////////////////////////////////////////////////////
    // STAGE 1: Arrange arc storage for each node.
    for (p_node = &*m_nodes.begin(); p_node != &*m_nodes.end(); ++p_node) {

        size_type fwd_arc_count
                = reinterpret_cast<size_type>(p_node->p_first_out_arc);
        size_type rev_arc_count
                = reinterpret_cast<size_type>(p_node->p_first_in_arc );

        assert(fwd_arc_count < m_fwd_arcs.chunk_size());
        assert(rev_arc_count < m_rev_arcs.chunk_size());

        // Process forward arcs.
        if (p_fwd_arc + fwd_arc_count > p_fwd_chunk->end()) {

            // If the current chunk is full, proceed to the next chunk.
            if (fwd_is_full
            || (++it_fwd_chunk_ptr == m_fwd_arcs.list().end())) {
                // Begin allocation of new chunks.
                p_fwd_chunk   = m_fwd_arcs.new_chunk();
                fwd_is_full   = true;
                p_rev_arc1    = 0;
            }
            else {
                p_fwd_chunk   = *it_fwd_chunk_ptr;
                p_rev_arc1    = (*(++it_rev_chunk_ptr1))->begin();
            }

            p_fwd_arc         = p_fwd_chunk->begin();
            p_node1           = p_node - 1; // The previous node.
            p_node1->tag     |= SPECIAL_OUT;
            p_node1->dist_id  = reinterpret_cast<size_t>
                                    (reinterpret_cast<fwd_arc*>
                                        (p_node1->dist_id) + 1);

            if (p_node1->p_parent_arc) {
                p_node1->p_parent_arc = reinterpret_cast<fwd_arc*>
                                            (reinterpret_cast<rev_arc*>
                                                (p_node1->p_parent_arc) + 1);
            }

            ++(p_node1->p_first_out_arc);
        }

        // Process reverse arcs.
        if (p_rev_arc + rev_arc_count > p_rev_chunk->end()) {

            // If the current chunk is full, proceed to the next chunk.
            if (rev_is_full 
            || (++it_rev_chunk_ptr == m_rev_arcs.list().end())) {
                // Begin allocation of new chunks.
                p_rev_chunk = m_rev_arcs.new_chunk();
                rev_is_full = true;
            }
            else {
                p_rev_chunk = *it_rev_chunk_ptr;
            }

            p_rev_arc       = p_rev_chunk->begin();
            p_node1         = p_node - 1; // The previous node.
            p_node1->tag   |= SPECIAL_IN;
            p_node1->dist   = reinterpret_cast<size_t>
                                  (reinterpret_cast<rev_arc*>
                                      (p_node1->dist) + 1);

            ++(p_node1->p_first_in_arc);
        }

        // Increment the arc pointers by the arc counts of the current node.
        p_fwd_arc += fwd_arc_count;
        p_rev_arc += rev_arc_count;

        if (p_rev_arc1 != 0) {
            p_rev_arc1           += fwd_arc_count;
            p_node->p_parent_arc  = reinterpret_cast<fwd_arc*>(p_rev_arc1);
        }
        else {
            p_node->p_parent_arc  = 0;
        }

        // Points to the last outgoing arc of the node.
        p_node->p_first_out_arc   = p_fwd_arc;
        // Points to the last incoming arc of the node.
        p_node->p_first_in_arc    = p_rev_arc;

        p_node->dist_id           = reinterpret_cast<size_t>(p_fwd_arc);
        p_node->dist              = reinterpret_cast<size_t>(p_rev_arc);

    } // for

    
    // The last node.
    p_node1                       = p_node - 1;
    p_node1->tag                 |= (SPECIAL_OUT | SPECIAL_IN);
    p_node1->dist_id              = reinterpret_cast<size_t>
                                        (reinterpret_cast<fwd_arc*>
                                            (p_node1->dist_id) + 1);
    p_node1->dist                 = reinterpret_cast<size_t>
                                        (reinterpret_cast<rev_arc*>
                                            (p_node1->dist) + 1);

    if (p_node1->p_parent_arc) {
        p_node1->p_parent_arc     = reinterpret_cast<fwd_arc*>
                                        (reinterpret_cast<rev_arc*>
                                            (p_node1->p_parent_arc) + 1);
    }

    ++(p_node1->p_first_out_arc);
    ++(p_node1->p_first_in_arc );

    ///////////////////////////////////////////////////////////////////////
    // STAGE 2: Sort forward arcs.
    //
    // Pre-condition:
    //     1) The arc slots are re-dispatched to every node (STEP 1).
    //     2) p_rev_arc->p_fwd = <pointer to tail node>
    //     3) p_fwd_arc->shift = <pointer to head node>
    //
    // Post-condition:
    //     1) The forward arcs are sorted in the order of their tail
    //        nodes.
    //
    for (it_fwd_chunk_ptr  = m_fwd_arcs.list().begin(),
         it_rev_chunk_ptr  = m_rev_arcs.list().begin(); 
         it_fwd_chunk_ptr != it_fwd_chunk_ptr_end;
         ++it_fwd_chunk_ptr,
         ++it_rev_chunk_ptr) {

        p_fwd_chunk = *it_fwd_chunk_ptr;
        p_rev_chunk = *it_rev_chunk_ptr;

        for (p_fwd_arc  = p_fwd_chunk->first(),
             p_rev_arc  = p_rev_chunk->first(); 
             p_fwd_arc != p_fwd_chunk->last();
             ++p_fwd_arc, ++p_rev_arc) {

            assert(p_fwd_arc != 0);
            assert(p_rev_arc != 0);

            capacity_type rev_cap_save, rev_cap = 0;
            ptrdiff_t     shift_save, shift = 0;

            // Get the tail node of the arc. It is NULL when the arc
            // is already processed, then we skip it.
            p_tail_node = reinterpret_cast<node*>(p_rev_arc->p_fwd);
            if (0 == p_tail_node)
                continue;

            p_rev_arc1 = p_rev_arc;
            p_fwd_arc1 = p_fwd_arc;

            do {

                p_rev_arc1->p_fwd = 0;

                // [shift] = [head node] - [tail node]
                shift_save   = p_fwd_arc1->shift - 
                                   reinterpret_cast<ptrdiff_t>(p_tail_node);
                // Save reverse residual capacity.
                rev_cap_save = p_fwd_arc1->rev_cap;

                if (shift != 0) {
                    p_fwd_arc1->rev_cap = rev_cap;
                    p_fwd_arc1->shift   = shift;
                }

                shift      = shift_save;
                rev_cap    = rev_cap_save;
                p_fwd_arc1 = --(p_tail_node->p_first_out_arc); // [!]

                // Move the arc information of the tail node to a new slot
                // allocated to the node. This new slot may be occupied as
                // well. We find out the tail node of the arc occuping the
                // new slot, and move the arc to the new slot allocated to
                // its tail node.  Keep this process until we reach an un-
                // occupied new slot.
                if (p_tail_node->p_parent_arc == 0) break;
                else {
                    p_tail_node->p_parent_arc 
                                = reinterpret_cast<fwd_arc*>
                                      (reinterpret_cast<rev_arc*>
                                          (p_tail_node->p_parent_arc) - 1);
                    p_rev_arc1  = reinterpret_cast<rev_arc*>
                                      (p_tail_node->p_parent_arc);
                }

            } while (0 != (
                            p_tail_node = reinterpret_cast<node *>
                                (p_rev_arc1->p_fwd)
                        )
                    );

            p_fwd_arc1->rev_cap = rev_cap;
            p_fwd_arc1->shift   = shift;
            
        } // for
    }

    ///////////////////////////////////////////////////////////////////////
    // STAGE 3: Sort reverse arcs.
    for (p_node = &*m_nodes.begin(); p_node != &*m_nodes.end(); ++p_node) {
        for (p_fwd_arc = p_node->p_first_out_arc; 
             p_fwd_arc < reinterpret_cast<fwd_arc*>(p_node->dist_id); 
             ++p_fwd_arc) {

            p_head_node = reinterpret_cast<node*>
                              (reinterpret_cast<char*>
                                  (p_node) + p_fwd_arc->shift);
            (--p_head_node->p_first_in_arc)->p_fwd = p_fwd_arc;
        }
    }

    ///////////////////////////////////////////////////////////////////////
    // STAGE 4: Process special nodes.
    for (p_node = &*m_nodes.begin(); p_node != &*m_nodes.end(); ++p_node) {
        if (p_node->tag & SPECIAL_OUT) {
            // Special out arc : shift = p_last_out_arc
            (--p_node->p_first_out_arc)->shift = p_node->dist_id;
        }
        if (p_node->tag & SPECIAL_IN) {
            // Special in arc  : p_fwd = p_last_in_arc
            (--p_node->p_first_in_arc)->p_fwd 
                = reinterpret_cast<fwd_arc*>(p_node->dist);
        }
    }
    
    m_prepared = true;
}

//
// Constants
//
template<typename _Cap, typename _Tg>
    const unsigned char graph_fs<_Cap, _Tg>::SPECIAL_IN  = 0x20;

template<typename _Cap, typename _Tg>
    const unsigned char graph_fs<_Cap, _Tg>::SPECIAL_OUT = 0x40;

} // namespace

#endif
