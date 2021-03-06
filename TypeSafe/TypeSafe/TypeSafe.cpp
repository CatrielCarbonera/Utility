// TypeSafe.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <string>

#include "Utility.h"
DECLARE_SAFE_TYPE(AngleType, double);
DECLARE_SAFE_TYPE(FileNameType, std::string);
DECLARE_SAFE_TYPE(PersonNameType, std::string);

void print(AngleType const & a)
{
    std::cout << a.get() << std::endl;
}

void print(FileNameType const & a)
{
    std::cout << a.get() << std::endl;
}

void print(PersonNameType const & a)
{
    std::cout << a.get() << std::endl;
}

int main()
{
    AngleType a(2.);
    AngleType b(3);
    print(a + b);
    print(a - b);
    print(a * 3.5);
    
    FileNameType fn(__FILE__);
    print(fn);

    PersonNameType my_name("Catriel");

#ifdef COMPILE_ERROR
    print(PersonNameType(fn));
#endif
    char const * str1 = "string type 1";
    char const * str2 = "string type 2";

    return 0;
}

