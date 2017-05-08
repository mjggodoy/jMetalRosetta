// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_sicdock_xyzStripeHashWithMeta_fwd_hh
#define INCLUDED_protocols_sicdock_xyzStripeHashWithMeta_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <numeric/types.hh>

namespace numeric {
namespace geometry {
namespace hashing {

template<typename T>
class xyzStripeHashWithMeta;
typedef utility::pointer::shared_ptr< xyzStripeHashWithMeta<numeric::Real> > xyzStripeHashWithMetaRealOP;
typedef utility::pointer::shared_ptr< xyzStripeHashWithMeta<numeric::Real> const > xyzStripeHashWithMetaRealCOP;
typedef utility::pointer::shared_ptr< xyzStripeHashWithMeta<float> > xyzStripeHashWithMetaFloatOP;
typedef utility::pointer::shared_ptr< xyzStripeHashWithMeta<float> const > xyzStripeHashWithMetaFloatCOP;

}
}
}

#endif
