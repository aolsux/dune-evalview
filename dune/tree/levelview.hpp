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

#pragma once
#include <tree/root.hpp>


namespace tree {

template< class GV > class Root;
    
template< class GV > 
class LevelView {
public:
  typedef foo Node;
  std::vector<Node*> _nodes;
  
  
  // create an iterator that resolves the double dereferencation to a single one
  struct const_iterator : public std::vector<Node*>::const_iterator
  {
      typedef std::vector<Node*>::iterator Base;
      
      iterator(std::vector<Node*>::iterator it) : Base(it)
      {}
      
      const Node& operator *() const
      {
        return *Base::operator*();
      }
  };
  
  // make the view const iterable
  const_iterator begin() const     { return const_iterator(_nodes.begin());}
  const_iterator end() const       { return const_iterator(_nodes.end());}
  
private:
    LevelView() = delete;
protected:
  
   // NOTE: the level argument is a possible candidate for an additional template argument
    LevelView(const Root<GV>& root, unsigned level)
    {
       recurse(root,level);
    }
    
    
    void recurse(const Node& node, unsigned level)
    {
       if(node.level() == level)
         _nodes.push_back(&node);
       else
         //NOTE: node.child(0) returns a pointer, or NULL if no child exists
         if(node.child(0) != NULL) recurse(*node.child(0));
         if(node.child(1) != NULL) recurse(*node.child(1));
       assert(node.level()<=level,"improper recursion limit criterion");
    }
    
    
      
    
};


}
