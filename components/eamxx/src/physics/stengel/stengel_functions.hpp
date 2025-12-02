#ifndef STENGEL_FUNCTIONS_HPP
#define STENGEL_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace stengel {

template <typename ScalarT, typename DeviceT>
struct StengelFunctions
{

  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using Mask = ekat::Mask<BigPack<Scalar>::n>;
  using Smask = ekat::Mask<SmallPack<Scalar>::n>;

  // GPU/Kokkos related things
  using KT = ekat::KokkosTypes<Device>;
  using TeamPolicy  = typename KokkosTypes<Device>::TeamPolicy;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  // KOKKOS_FUNCTION 
  static void scale_field_2d(view_2d<Spack> &field, Real alpha);
  // KOKKOS_FUNCTION 
  static void get_max_value(view_2d<Spack> &field, Real &max_value);

}; // struct Functions

template<typename S, typename D>
// KOKKOS_FUNCTION
void StengelFunctions<S,D>::scale_field_2d(view_2d<Spack> &field, Real alpha) {

  const int ni = static_cast<int>(field.extent(0));
  const int nj = static_cast<int>(field.extent(1));

  using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
  MDPolicy mdp({0,0}, {ni,nj});

  // create a device mirror and copy host->device, run kernel, copy device->host
  auto field_dev = Kokkos::create_mirror_view_and_copy(DefaultDevice(), field); // or try HostDevice()
  Kokkos::parallel_for(mdp, KOKKOS_LAMBDA(const int i, const int j) {
      field_dev(i,j) *= alpha;
  });
  Kokkos::deep_copy(field, field_dev);
};

template<typename S, typename D>
// KOKKOS_FUNCTION
void StengelFunctions<S,D>::get_max_value(view_2d<Spack> &field, Real &max_value) {
  // compute global max over all i,j in a single parallel_reduce (no nested lambdas)
  const int ni = static_cast<int>(field.extent(0));
  const int nj = static_cast<int>(field.extent(1));
  Real global_max = -std::numeric_limits<Real>::infinity();

  using MDPolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
  MDPolicy mdp({0,0}, {ni,nj});

  Kokkos::parallel_reduce("GetMaxValueMD", mdp,
    KOKKOS_LAMBDA(const int i, const int j, Real &local_max) {
      const Real v = field(i,j)[0];
      if (v > local_max) local_max = v;
    },
    Kokkos::Max<Real>(global_max)
  );

  // store result on host
  max_value = global_max;
};

} // namespace stengel
} // namespace scream

#endif // STENGEL_FUNCTIONS_HPP
