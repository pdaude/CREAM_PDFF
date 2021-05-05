/*
 ==========================================================================
 |   
 |   $Id: chunk_list.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___CHUNK_LIST_HXX___
#   define ___CHUNK_LIST_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4571)
#   endif

#   include <cstddef>   // for size_t
#   include <cstring>   // for ::memset
#   include <list>      // for std::list

#   include <optnet/_base/iterator.hxx>

///  @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class chunk_list
///  @brief A list that allocates storage in chunks.
///////////////////////////////////////////////////////////////////////////
template <class _Ty, size_t _ChunkSize = 1024>
class chunk_list
{
    typedef chunk_list<_Ty, _ChunkSize> _Self;

public:

    ///////////////////////////////////////////////////////////////////////
    ///  @class chunk
    ///  @brief The class that represent a chunk of elements of the list.
    ///////////////////////////////////////////////////////////////////////
    struct chunk
    {
        chunk()
        {
            ::memset(m_adata, 0, sizeof(m_adata));
            m_plast = &(m_adata[0]);
        }

        inline _Ty* begin() { return &(m_adata[0]); }
        inline _Ty* end()   { return &(m_adata[_ChunkSize]); }
        inline _Ty* first() { return &(m_adata[0]); }
        inline _Ty* last()  { return m_plast; }

        _Ty* m_plast;
        _Ty  m_adata[_ChunkSize + 1];
    };

    typedef std::list<chunk*>   list_type;
    typedef chunk*              chunk_pointer;
    typedef const chunk*        chunk_const_pointer;

    typedef _Ty                 value_type;
    typedef value_type*         pointer;
    typedef const value_type*   const_pointer;
#   ifndef __OPTNET_CRAPPY_MSC__
    typedef optnet::detail::normal_iterator<pointer, _Self>
                                iterator;
    typedef optnet::detail::normal_iterator<const_pointer, _Self>
                                const_iterator;
#   else
    typedef chunk*              iterator;
    typedef const chunk*        const_iterator;
#   endif
    typedef value_type&         reference;
    typedef const value_type&   const_reference;
    typedef size_t              size_type;


    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///////////////////////////////////////////////////////////////////////
    chunk_list() : m_n(0)
    {
        // Do nothing.
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Default destructor.
    ///////////////////////////////////////////////////////////////////////
    ~chunk_list()
    {
        clear();
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Append an element to the end of the list.
    ///
    ///  @return A pointer to the newly-appended element.
    ///
    ///  @remarks This method allocates a new chunk if needed.
    ///
    ///  @see new_chunk
    ///
    ///////////////////////////////////////////////////////////////////////
    pointer append()
    {
        pointer p;

        if (m_list.empty() || 
            m_list.back()->last() >= m_list.back()->end()) {
            
            // Reached the end of a chunk, so create a new one.
            new_chunk();

            // Return pointer to the new slot.
            p = m_list.back()->m_plast++;
            ++m_n;

        } else {

            // Return pointer to the new slot.
            p = m_list.back()->m_plast++;
            ++m_n;
        }

        return p;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Allocate a new chunk and append it to the end of the chunk list.
    ///
    ///  @returns A pointer to the newly-created chunk object.
    ///////////////////////////////////////////////////////////////////////
    chunk_pointer new_chunk()
    {
        chunk_pointer pchunk = new chunk();
        m_list.push_back(pchunk);
        return pchunk;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Frees all elements in the list.
    ///////////////////////////////////////////////////////////////////////
    void clear()
    {
        typename list_type::iterator it;
        for (it = m_list.begin(); it != m_list.end(); ++it)
            delete *it;
        m_list.clear();
        m_n = 0;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the number of elements in the list.
    ///////////////////////////////////////////////////////////////////////
    inline size_type    count() const           { return m_n; }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the size of a chunk, i.e., the number of elements that
    ///  a chunk can accomodate.
    ///////////////////////////////////////////////////////////////////////
    inline size_type    chunk_size() const      { return _ChunkSize; }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a list of type std::list<chunk*> that contains all chunks
    ///  in the list.
    ///////////////////////////////////////////////////////////////////////
    inline list_type&   list()                  { return m_list; }


private:

    list_type m_list;
    size_type m_n;
};


} // namespace

#endif
