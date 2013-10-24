
#pragma once

#include <root.hpp>


namespace evalview {

class LeafView {
  LeafView() = delete;
protected:
  LeafView(const Root& root);
};
}
