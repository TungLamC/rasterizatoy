#include <variant>
#include <iostream>

#include "rasterizatoy.hpp"

using namespace rasterizatoy;

typedef uintptr_t gsize;

int main() {
	Matrix4D matrix {
		1, 2, 3, 4,
		1, 5, 3, 4,
		1, 2, 7, 4,
		1, 2, 3, 4,
	};
	std::cout << matrix << std::endl;
	return 0;
}

