#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

int main() {
    namespace fs = std::filesystem;
    const fs::path here = fs::current_path();
    const fs::path exe = here / "litio_shooting";
    if (!fs::exists(exe)) {
        std::cerr << "No existe ejecutable: " << exe << "\n";
        return 1;
    }

    int rc = std::system("./litio_shooting");
    if (rc != 0) {
        return rc;
    }

    // Guarda copia con etiqueta C++ para trazabilidad.
    const fs::path src_e = here / "energias_shooting_fortran.dat";
    const fs::path src_w = here / "funciones_shooting_fortran.dat";
    if (fs::exists(src_e)) {
        fs::copy_file(src_e, here / "energias_shooting_cpp.dat", fs::copy_options::overwrite_existing);
    }
    if (fs::exists(src_w)) {
        fs::copy_file(src_w, here / "funciones_shooting_cpp.dat", fs::copy_options::overwrite_existing);
    }
    return 0;
}
