module ModUser

  use ModUserEmpty

  include 'user_module.h'

  character(len=*), parameter :: NameUserFile = 'srcUser/ModUserEarlyEarthSkeleton.f90'
  character(len=*), parameter :: NameUserModule = 'EARLY EARTH SKELETON (EMPTY HOOKS)'

end module ModUser
