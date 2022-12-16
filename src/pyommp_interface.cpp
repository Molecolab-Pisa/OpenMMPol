#include "openmmpol.h"
#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;


class OMMPSystem{
    public:
        OMMPSystem(std::string mmp_filename) : handler(nullptr){
            handler = ommp_init_mmp(mmp_filename.c_str());
        }

        ~OMMPSystem(){
            if(is_init()) ommp_terminate(handler);
            handler = nullptr;
        }
        
        bool is_init(){
            return handler != nullptr;
        }

        void print_summary(){
            ommp_print_summary(handler);
        }

    private: 
        void *handler;
};

PYBIND11_MODULE(pyopenmmpol, m){
    py::class_<OMMPSystem, std::shared_ptr<OMMPSystem>>(m, "OMMPSystem", "System of OMMP library.")
        .def(py::init<std::string>(), "provola", "mmp_filename")
        .def("print_summary", &OMMPSystem::print_summary, "provolotto");
}
