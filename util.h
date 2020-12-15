#ifndef __util_h__
#define __util_h__

#include "parser.h"

namespace util
{

    class Util
    {
    private:
        friend std::ostream& operator<<(std::ostream& os, const parser::Vec3f& vec);
        friend std::ostream& operator<<(std::ostream& os, const parser::Vec3i& vec);
        friend std::ostream& operator<<(std::ostream& os, const parser::Vec4f& plane);
        friend std::ostream& operator<<(std::ostream& os, const parser::Camera& cam);
        friend std::ostream& operator<<(std::ostream& os, const parser::PointLight& pl);
        friend std::ostream& operator<<(std::ostream& os, const parser::Material& mat);
        friend std::ostream& operator<<(std::ostream& os, const parser::Sphere& sph);
        friend std::ostream& operator<<(std::ostream& os, const parser::Scene& scene);

    public:
        void PrintSceneDetails(parser::Scene scene);

    };

}
#endif // __util_h__