#ifndef COMMON_H
#define COMMON_H

#include "proj.hpp"
using namespace Projection;


typedef struct {
    const char* name;
    enum Projection proj;
    const char* params;
    Center center;
    Layout layout;
    SquareMethod square_method;
} projection_t;

projection_t projections[] = {
#if 1
    // Hemisphere Layouts
    { "lea-hemisphere-north",    ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    //{ "lea-hemisphere-south",    ProjLambertAzimuthalEqualArea,  "", CenterSouthPole, LayoutHemisphere, SquareMethodNone },
    { "stg-hemisphere-north",     ProjStereographic,              "", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    //{ "stg-hemisphere-south",     ProjStereographic,              "", CenterSouthPole, LayoutHemisphere, SquareMethodNone },
    { "eqd-hemisphere-north",     ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    //{ "eqd-hemisphere-south",     ProjEquidistantAzimuthal,       "", CenterSouthPole, LayoutHemisphere, SquareMethodNone },
    { "bhm-hemisphere-north",     ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "bhm-hemisphere-south",     ProjBreusingHarmonicMean,       "", CenterSouthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    //{ "mix-hemisphere-south",     ProjMixLambertStereographic,    "beta=0.4", CenterSouthPole, LayoutHemisphere, SquareMethodNone },
#endif
#if 1
    // Squared Hemisphere Layouts
    //{ "lea-hemisphere-shirley-north", ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutHemisphere, SquareMethodShirley },
    //{ "lea-hemisphere-shirley-south", ProjLambertAzimuthalEqualArea,  "", CenterSouthPole, LayoutHemisphere, SquareMethodShirley },
    //{ "stg-hemisphere-conf-north",     ProjStereographic,              "", CenterNorthPole, LayoutHemisphere, SquareMethodConformal },
    //{ "stg-hemisphere-conf-south",     ProjStereographic,              "", CenterSouthPole, LayoutHemisphere, SquareMethodConformal },
    //{ "eqd-hemisphere-ell-north",      ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutHemisphere, SquareMethodElliptical },
    //{ "eqd-hemisphere-ell-south",      ProjEquidistantAzimuthal,       "", CenterSouthPole, LayoutHemisphere, SquareMethodElliptical },
    { "bhm-hemisphere-ell-north",      ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutHemisphere, SquareMethodElliptical },
    { "bhm-hemisphere-ell-south",      ProjBreusingHarmonicMean,       "", CenterSouthPole, LayoutHemisphere, SquareMethodElliptical },
    //{ "mix-hemisphere-sq-north",       ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutHemisphere, SquareMethodSquircle },
    //{ "mix-hemisphere-sq-south",       ProjMixLambertStereographic,    "beta=0.4", CenterSouthPole, LayoutHemisphere, SquareMethodSquircle },
#endif
#if 1
    // Sphere Layouts
    { "lea-sphere-north",    ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutSphere, SquareMethodNone },
    //{ "lea-sphere-south",    ProjLambertAzimuthalEqualArea,  "", CenterSouthPole, LayoutSphere, SquareMethodNone },
    { "eqd-sphere-north",     ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutSphere, SquareMethodNone },
    //{ "eqd-sphere-south",     ProjEquidistantAzimuthal,       "", CenterSouthPole, LayoutSphere, SquareMethodNone },
    { "bhm-sphere-north",     ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutSphere, SquareMethodNone },
    //{ "bhm-sphere-south",     ProjBreusingHarmonicMean,       "", CenterSouthPole, LayoutSphere, SquareMethodNone },
#endif
#if 1
    // Squared Sphere Layouts
    //{ "lea-sphere-shirley-north", ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutSphere, SquareMethodShirley  },
    //{ "lea-sphere-shirley-south", ProjLambertAzimuthalEqualArea,  "", CenterSouthPole, LayoutSphere, SquareMethodShirley  },
    { "eqd-sphere-sq-north",       ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutSphere, SquareMethodSquircle },
    //{ "eqd-sphere-sq-south",       ProjEquidistantAzimuthal,       "", CenterSouthPole, LayoutSphere, SquareMethodSquircle },
    //{ "bhm-sphere-sq-north",       ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutSphere, SquareMethodSquircle },
    //{ "bhm-sphere-sq-south",       ProjBreusingHarmonicMean,       "", CenterSouthPole, LayoutSphere, SquareMethodSquircle },
#endif
#if 1
    // Quincuncial Layouts
    //{ "lea-quinc-stretch",     ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutQuincuncial, SquareMethodStretch    },
    //{ "lea-quinc-adjstretch",  ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutQuincuncial, SquareMethodAdjStretch },
    { "lea-quinc-shirley",     ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutQuincuncial, SquareMethodShirley    },
    //{ "lea-quinc-squircle",    ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutQuincuncial, SquareMethodSquircle   },
    //{ "lea-quinc-elliptic",    ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutQuincuncial, SquareMethodElliptical },
    //{ "lea-quinc-conf",        ProjLambertAzimuthalEqualArea,  "", CenterNorthPole, LayoutQuincuncial, SquareMethodConformal  },
    //{ "stg-quinc-stretch",      ProjStereographic,              "", CenterNorthPole, LayoutQuincuncial, SquareMethodStretch    },
    //{ "stg-quinc-adjstretch",   ProjStereographic,              "", CenterNorthPole, LayoutQuincuncial, SquareMethodAdjStretch },
    //{ "stg-quinc-shirley",      ProjStereographic,              "", CenterNorthPole, LayoutQuincuncial, SquareMethodShirley    },
    //{ "stg-quinc-squircle",     ProjStereographic,              "", CenterNorthPole, LayoutQuincuncial, SquareMethodSquircle   },
    //{ "stg-quinc-elliptic",     ProjStereographic,              "", CenterNorthPole, LayoutQuincuncial, SquareMethodElliptical },
    { "stg-quinc-conf",         ProjStereographic,              "", CenterNorthPole, LayoutQuincuncial, SquareMethodConformal  },
    //{ "mix-quinc-stretch",      ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutQuincuncial, SquareMethodStretch    },
    //{ "mix-quinc-adjstretch",   ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutQuincuncial, SquareMethodAdjStretch },
    //{ "mix-quinc-shirley",      ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutQuincuncial, SquareMethodShirley    },
    { "mix-quinc-squircle",     ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutQuincuncial, SquareMethodSquircle   },
    //{ "mix-quinc-elliptic",     ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutQuincuncial, SquareMethodElliptical },
    //{ "mix-quinc-conf",         ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutQuincuncial, SquareMethodConformal  },
    //{ "bhm-quinc-stretch",      ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodStretch    },
    //{ "bhm-quinc-adjstretch",   ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodAdjStretch },
    //{ "bhm-quinc-shirley",      ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodShirley    },
    //{ "bhm-quinc-squircle",     ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodSquircle   },
    { "bhm-quinc-elliptic",     ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodElliptical },
    //{ "bhm-quinc-conf",         ProjBreusingHarmonicMean,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodConformal  },
    //{ "eqd-quinc-stretch",      ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodStretch    },
    //{ "eqd-quinc-adjstretch",   ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodAdjStretch },
    //{ "eqd-quinc-shirley",      ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodShirley    },
    //{ "eqd-quinc-squircle",     ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodSquircle   },
    { "eqd-quinc-elliptic",     ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodElliptical },
    //{ "eqd-quinc-conf",         ProjEquidistantAzimuthal,       "", CenterNorthPole, LayoutQuincuncial, SquareMethodConformal  },
#endif
#if 1
    { "mix-00-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.0", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-01-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.1", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-02-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.2", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-03-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.3", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-04-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.4", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-05-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.5", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-06-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.6", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-07-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.7", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-08-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.8", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-09-hemisphere-north",     ProjMixLambertStereographic,    "beta=0.9", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
    { "mix-10-hemisphere-north",     ProjMixLambertStereographic,    "beta=1.0", CenterNorthPole, LayoutHemisphere, SquareMethodNone },
#endif
    { "perspective", ProjPerspective, "", CenterNorthPole, LayoutHemisphere, SquareMethodNone }
};

int n_projections = sizeof(projections) / sizeof(projections[0]);

#endif
