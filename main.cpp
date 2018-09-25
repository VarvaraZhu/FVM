#include<iostream>
#pragma once

int main(int argc, char* argv[])
{
       if (argc > 1) std::cout << "Mesh will be read from " << argv[1] << std::endl;
       else {
                 std::cout << "Please, set the mesh file" << std::endl;
                 return 0;
       }




       return 0;
}
