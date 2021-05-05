/*
 ==========================================================================
 |   
 |   $Id: optnet_fs_3d_multi.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___OPTNET_FS_3D_MULTI_HXX___
#   define ___OPTNET_FS_3D_MULTI_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_fs/optnet_fs_maxflow.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <vector>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_fs_3d_multi
///  @brief Implementation of the multiple-surface Optimal Net algorithm
///         using the Boykov-Kolmogorov max-flow algorithm on a
///         forward-star represented graph.
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg = net_f_xy>
class optnet_fs_3d_multi
{
    typedef optnet_fs_maxflow<_Cap, net_f_xy>   graph_type;

    struct  _Intra {

        std::vector<                    //
            std::pair<size_t, int>      // A tiny adjacency graph
        >       tonode[2];              //
        bool    circle[2];
        int     margin[2];
        int     smooth[2];
        int     visits[2];
        bool    dropped;
        bool    climbed;
    };

    struct  _Inter {
        size_t  k[2];
        int     r[2];
    };

    typedef std::vector<_Intra>                 intra_vector;
    typedef std::vector<_Inter>                 inter_vector;


public:

    typedef typename graph_type::roi_base_type  roi_base_type;
    typedef typename graph_type::roi_ref_type   roi_ref_type;
    typedef typename graph_type::roi_type       roi_type;

    typedef size_t                              size_type;
    typedef _Cost                               cost_type;
    typedef _Cap                                capacity_type;

    typedef array_base<cost_type, _Tg>          cost_array_base_type;
    typedef array_ref<cost_type, _Tg>           cost_array_ref_type;
    typedef array<cost_type, _Tg>               cost_array_type;
    
    typedef array_base<size_type>               net_base_type;
    typedef array_ref<size_type>                net_ref_type;
    typedef array<size_type>                    net_type;
    
    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    optnet_fs_3d_multi();

    ///////////////////////////////////////////////////////////////////////
    ///  Create graph and set the cost values.
    ///
    ///  @param cost The array that contain the cost values for each voxel.
    ///
    ///  @remarks This function will re-assign the node costs of the
    ///           underlying graph based upon the given cost array.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create(const cost_array_base_type& cost);

    ///////////////////////////////////////////////////////////////////////
    ///  Set the smoothness constraints for the k-th surface.
    ///
    ///  @param k       The index to the k-th surface.
    ///  @param smooth0 The first smoothness parameter.
    ///  @param smooth1 The second smoothness parameter.
    ///  @param circle0 Enabling/disabling circle graph construction.
    ///  @param circle1 Enabling/disabling circle graph construction.
    ///
    ///  @remarks The actual meanings of the smoothness parameters depend
    ///           on the direction setting of the net. 
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_params(size_type   k,
                    int         smooth0 = 1,
                    int         smooth1 = 1,
                    bool        circle0 = false,
                    bool        circle1 = false
                    );


    ///////////////////////////////////////////////////////////////////////
    ///  Set the relation between a pair of surfaces.
    ///
    ///  @param  k0  The index to one surface.
    ///  @param  k1  The index to the other surface.
    ///  @param  r01 The relation between the two surfaces.
    ///  @param  r10 The relation between the two surfaces.
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_relation(size_type k0, size_type k1, int r01, int r10);

#   ifdef __OPTNET_SUPPORT_ROI__

    ///////////////////////////////////////////////////////////////////////
    ///  Set the region of interest.
    ///
    ///  @param roi The region-of-interest mask array.
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_roi(const roi_base_type& roi);

#   endif // __OPTNET_SUPPORT_ROI__

    ///////////////////////////////////////////////////////////////////////
    ///  Solve the optimal surface problem using the given cost function
    ///  and smoothness constraints.
    ///
    ///  @param net   The resulting optimal "net" surface.
    ///  @param pflow The output maximum flow value.
    ///
    ///  @remarks The pflow parameter, if not NULL, will return the
    ///           computed maximum-flow value. It is used primarily
    ///           for debugging.
    ///
    ///////////////////////////////////////////////////////////////////////
    void solve(net_base_type& net,      // [OUT]
               capacity_type* pflow = 0 // [OUT]
               );


private:

    ///////////////////////////////////////////////////////////////////////
    // Compute the upper and lower margin of the 3-D subgraphs. The nodes
    // above the upper bound and below the lower bound can never be on
    // the resulting surfaces.
    void get_bounds_of_subgraphs();

    ///////////////////////////////////////////////////////////////////////
    // Transform the costs of the graph nodes based on the given
    // cost vector.
    void transform_costs();
    
    ///////////////////////////////////////////////////////////////////////
    // Construct the arcs of the underlying graph.
    void build_arcs();

    ///////////////////////////////////////////////////////////////////////
    const cost_array_base_type* m_pcost;

    graph_type      m_graph;
    intra_vector    m_intra;
    inter_vector    m_inter;
};

} // optnet

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_fs/optnet_fs_3d_multi.cxx>
#   endif

#endif
