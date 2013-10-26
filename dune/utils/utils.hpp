//**************************************************************************************//
//     AUTHOR: Malik Kirchner "malik.kirchner@gmx.net"                                  //
//             Martin Rückl "martin.rueckl@physik.hu-berlin.de"                         //
//                                                                                      //
//     This program is free software: you can redistribute it and/or modify             //
//     it under the terms of the GNU General Public License as published by             //
//     the Free Software Foundation, either version 3 of the License, or                //
//     (at your option) any later version.                                              //
//                                                                                      //
//     This program is distributed in the hope that it will be useful,                  //
//     but WITHOUT ANY WARRANTY; without even the implied warranty of                   //
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    //
//     GNU General Public License for more details.                                     //
//                                                                                      //
//     You should have received a copy of the GNU General Public License                //
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.            //
//                                                                                      //
//     Dieses Programm ist Freie Software: Sie können es unter den Bedingungen          //
//     der GNU General Public License, wie von der Free Software Foundation,            //
//     Version 3 der Lizenz oder (nach Ihrer Option) jeder späteren                     //
//     veröffentlichten Version, weiterverbreiten und/oder modifizieren.                //
//                                                                                      //
//     Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber           //
//     OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite               //
//     Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.       //
//     Siehe die GNU General Public License für weitere Details.                        //
//                                                                                      //
//     Sie sollten eine Kopie der GNU General Public License zusammen mit diesem        //
//     Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.       //
//                                                                                      //
//**************************************************************************************//
/*! \file */
#pragma once

#include <sstream>
#include <string>
#include <functional>
#include <sys/time.h>

// #include <error/baseerror.hpp>
// #include <utils/tuple.hpp>

template< typename _T >
inline void safe_delete( _T*& pnt ) {
    if ( pnt ) delete pnt;
    pnt = NULL;
}

template< typename _T >
inline void safe_array_delete( _T*& pnt ) {
    if ( pnt ) delete [] pnt;
    pnt = NULL;
}

template< typename T >
inline const std::string asString( const T& arg ) {
    std::stringstream buf;

    buf << arg;

    return buf.str();
}

inline const std::string asString( const std::string& arg ) {
    return arg;
}

inline const std::string asString( const bool arg ) {
    return arg ? "true" : "false";
}

enum CnslEsc {
    CE_RESET        = 0,    // Reset all attributes
    CE_BRIGHT       = 1,    // Bright
    CE_DIM          = 2,    // Dim
    CE_UNDERSCORE   = 4,    // Underscore
    CE_BLINK        = 5,    // Blink
    CE_REVERSE      = 7,    // Reverse
    CE_HIDDEN       = 8,    // Hidden

// Foreground Colors
    CE_FG_BLACK     = 30,   // Black
    CE_FG_RED       = 31,   // Red
    CE_FG_GREEN     = 32,   // Green
    CE_FG_YELLOW    = 33,   // Yellow
    CE_FG_BLUE      = 34,   // Blue
    CE_FG_MAGENTA   = 35,   // Magenta
    CE_FG_CYAN      = 36,   // Cyan
    CE_FG_WHITE     = 37,   // White

// Background Colors
    CE_BG_BLACK     = 40,   // Black
    CE_BG_RED       = 41,   // Red
    CE_BG_GREEN     = 42,   // Green
    CE_BG_YELLOW    = 43,   // Yellow
    CE_BG_BLUE      = 44,   // Blue
    CE_BG_MAGENTA   = 45,   // Magenta
    CE_BG_CYAN      = 46,   // Cyan
    CE_BG_WHITE     = 47    // White
};

inline std::ostream& operator << ( std::ostream& os, const CnslEsc& val ){
    char esc[10];
    sprintf(esc, "\033[%dm", val);
    os << esc;
    return os;
}

#define COLOR_CONSOLE
// enable colored console output with "-DCOLOR_CONSOLE"
#ifdef COLOR_CONSOLE
    #define CE_KEYWORD  CE_FG_BLUE << CE_BRIGHT
    #define CE_LINE     CE_FG_GREEN
    #define CE_STATUS   CE_FG_GREEN
    #define CE_TIME     CE_FG_YELLOW
    #define CE_DEBUG    CE_FG_MAGENTA << CE_BRIGHT << "[ DEBUG ]   " << CE_RESET << CE_FG_MAGENTA
    #define CE_ERROR    CE_FG_RED  << CE_BRIGHT << "[ ERROR ]   " << CE_RESET << CE_FG_RED
    #define CE_WARNING  CE_FG_CYAN << CE_BRIGHT << "[ WARNING ] " << CE_RESET << CE_FG_CYAN
    #define CE_LICENSE  CE_FG_YELLOW
#else
    #define CE_KEYWORD  ""
    #define CE_LINE     ""
    #define CE_STATUS   ""
    #define CE_TIME     ""
    #define CE_DEBUG    "[ DEBUG ]   "
    #define CE_ERROR    "[ ERROR ]   "
    #define CE_WARNING  "[ WARNING ] "
    #define CE_LICENSE  ""
#endif


class Timer {
protected:
    double ticTime;
    double tocTime;
public:
    Timer();

    void        tic();
    double      toc();
    std::string tocStr();
};

Timer::Timer(){
    tic();
}

void        Timer::tic(){
//    timeval t;
//    gettimeofday(&t,0);
    ticTime = (double)clock()/(double)CLOCKS_PER_SEC; //(double)t.tv_usec/1000.0;
}

double      Timer::toc(){
//    timeval t;
//    gettimeofday(&t,0);
    tocTime = (double)clock()/(double)CLOCKS_PER_SEC; //(double)t.tv_usec/1000.0;
    return tocTime - ticTime;
}
