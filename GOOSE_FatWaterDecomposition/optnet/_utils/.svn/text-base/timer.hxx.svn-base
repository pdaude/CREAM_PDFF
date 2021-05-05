/*
 ==========================================================================
 |   
 |   $Id: timer.hxx 2137 2007-07-02 03:26:31Z kangli $
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

#ifndef ___TIMER_HXX___
#   define ___TIMER_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#   endif

#   include <optnet/config.h>

#   if defined(__OPTNET_OS_WINNT__)
#       include <windows.h>
#   else
#       include <ctime>
#   endif


/// @namespace optnet
namespace optnet {
    
    /// @namespace optnet::utils
    /// @brief The namespace that contains utility classes and funtions.
    namespace utils {


#   if defined(__OPTNET_OS_WINNT__)

///////////////////////////////////////////////////////////////////////////
///  @class timer
///  @brief Simple timer class for measuring elapsed time.
///
///  Under Windows, this class uses the Win32 high-resolution timing APIs:
///  QueryPerformanceFrequency() and QueryPerformanceCounter(). Under the
///  other platforms, this class use the C Standard Library clock()
///  function.
///////////////////////////////////////////////////////////////////////////
class timer
{
public:

    ///////////////////////////////////////////////////////////////////////
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    timer()
    {
        LARGE_INTEGER cfreq;

        // Set the clock frequency if it is not yet set.
        QueryPerformanceFrequency(&cfreq);
        m_cfreq = (double)cfreq.QuadPart;
        
        // Get the current value of the high-res performance counter.
        QueryPerformanceCounter(&m_start);
    }
    
    ///////////////////////////////////////////////////////////////////////
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    void restart()
    {
        // Get the current value of the high-res performance counter.
        QueryPerformanceCounter(&m_start);
    }
  
    ///////////////////////////////////////////////////////////////////////
    ///  Returns elapsed time since the timer is created or restarted.
    ///
    ///  @return (double) Elapsed time in seconds.
    ///
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    double elapsed() const
    {
        LARGE_INTEGER   now;
        QueryPerformanceCounter(&now);
        return double(now.QuadPart - m_start.QuadPart) / m_cfreq;
    }


private:
    
    // Frequency setting is based on the hardware clock that does not
    // change between calling, so set this one only once.
    DOUBLE          m_cfreq;

    // The starting time stored as a LARGE_INTEGER.
    LARGE_INTEGER   m_start;

};

#   else // !__OPTNET_OS_WINNT__

///////////////////////////////////////////////////////////////////////////
///  @class timer
///  @brief Simple timer class for measuring elapsed time.
///
///  Under Windows, this class uses the Win32 high-resolution timing APIs:
///  QueryPerformanceFrequency() and QueryPerformanceCounter(). Under the
///  other platforms, this class use the C Standard Library clock()
///  function.
///////////////////////////////////////////////////////////////////////////
class timer
{
public:

    ///////////////////////////////////////////////////////////////////////
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    timer()             { m_start = clock(); }
    
    ///////////////////////////////////////////////////////////////////////
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    void restart()      { m_start = clock(); }
  
    ///////////////////////////////////////////////////////////////////////
    ///  Returns elapsed time since the timer is created or restarted.
    ///
    ///  @return (double) Elapsed time in seconds.
    ///
    ///  @post elapsed()==0
    ///////////////////////////////////////////////////////////////////////
    double elapsed() const
    {
        return double(clock() - m_start) / CLOCKS_PER_SEC;
    }

private:

    // The starting time stored as a clock_t.
    clock_t m_start;
};

#   endif // __OPTNET_OS_WINNT__

    } // namespace
} // namespace

#endif
