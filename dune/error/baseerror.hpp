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

#include <exception>
#include <string>

#define __ERROR_INFO__    __FUNCTION__, __FILE__, __LINE__

class BaseError : public std::exception {
protected:
    std::string func;
    std::string file;
    int         line;

    const std::string where() const {
        char s[20];
        sprintf(s,"%5.5d", line);
        return file + " at " + func + ":" + s;
    }

public:
    BaseError( const char* fc, const char* f, const int l ) noexcept : func(fc), file(f), line(l) {}

    virtual const char* what() const noexcept {
        std::string msg = "Error in " + where();
        return msg.c_str();
    }
};




class NotImplemented : public BaseError {
public:
    NotImplemented ( const char* fc, const char* f, const int l ) : BaseError( fc, f, l ) {}

    virtual const char* what() const noexcept {
        std::string msg = "The called function/method was not implemented yet! " + where();
        return msg.c_str();
    }
};
