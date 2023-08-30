using DrWatson
@quickactivate "Macro"


#########################################
############### POTENTIALS ##############
#########################################
function ComputePotential(r, type)
  if type == "MelMel"
      Rmm = 0.00124 # (repulsion in mm^2/day of melanophores from melanophores)
      rmm = 0.02 # (length scale in millimeters for repulsion of melanophores from melanophores — appears in exponential)

      return Rmm * exp.(-r ./ rmm).*(r .<= 0.2);
  elseif type == "MelXan"
      Rxm = 0.00274 # (repulsion of melanophores from xanthophores)
      rxm = 0.02 #

      return Rxm * exp.(-r ./ rxm);
  elseif type == "XanMel"
      Rmx = 0.00226 # (repulsion of xanthophores from melanophores)
      rmx = 0.02#
      Amx = 0.001956 # (attraction of xanthophores by melanophores)
      amx = 0.012#

      return Rmx * exp.(-r ./ rmx) - Amx * exp.(-r ./ amx);
  elseif type == "XanXan"
      Rxx = 0.00055 # (repulsion in mm^2/day of xanthophores from xanthophores)
      rxx = 0.011 #

      return Rxx * exp.(-r ./ rxx);
  else
      println("no potential specified")
  end
end

#########################################
##### EXTENDING MATRICES FOR FFT ########
#########################################
# For now this function only works for
# square matrices of dimension n x n
function extend(A, numrows)
  m = numrows
  n, ~ = size(A)
  At = Float64.(zeros(n + 2 * m, n + 2 * m))
  # set up centre
  At[m+1:n+m, m+1:n+m] = A # works

  # top row
  At[1:m, m+1:m+n] = A[n-m+1:n, 1:n]

  # bottom row
  At[m+n+1:m+n+m, m+1:m+n] = A[1:m, 1:n]

  # left column
  At[m+1:m+n, 1:m] = A[1:n, n-m+1:n]

  # right column
  At[m+1:m+n, m+n+1:m+n+m] = A[1:n, 1:m]

  # upper left corner
  At[1:m, 1:m] = A[n-m+1:n, n-m+1:n]

  # upper right corner
  At[1:m, m+n+1:m+n+m] = A[n-m+1:n, 1:m]

  # lower right corner
  At[m+n+1:m+n+m, m+n+1:m+n+m] = A[1:m, 1:m]

  # lower left corner
  At[m+n+1:m+n+m, 1:m] = A[1:m, n-m+1:n]

  return At
end



#########################################
############## SAVING DATA ##############
#########################################

function DensityToCSV(folder, name, density; counter = "")
  if counter != ""
    ctrstr = @sprintf "_%07d.csv" counter
    tmp = string(folder, name, ctrstr)
    writedlm(tmp, density)
  else
    writedlm(string(folder, name, ".csv"), density)
  end
end




#########################################
##### NUMERICAL SOLVERS O(Δx)=1,2,4 #####
#########################################
# numerical solvers
function EEuStep(t, y, p, Δt, RHS; parallel = false)
  k1 = RHS(t, y, p; parallel = parallel)
  y += Δt * k1
  t += Δt
  return (t, y)
end

function RK2Step(t, y, Δt, RHS; args...)
  k1 = RHS(t, y; args...)
  k2 = RHS(t + Δt / 2, y + k1 .* Δt / 2; args...)

  y += Δt * k2
  t += Δt
  return (t, y)
end

function RK4Step(t, y, Δt, RHS; args...)

  k1 = RHS(t, y; args...)
  k2 = RHS(t + Δt / 2, y + k1 .* Δt / 2; args...)
  k3 = RHS(t + Δt / 2, y + k2 .* Δt / 2; args...)
  k4 = RHS(t + Δt, y + k3 .* Δt; args...)

  y += Δt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
  t += Δt
  return (t, y)
end
#########################################
#########################################
#########################################





#########################################
####### SETTING UP DOMAIN & MESH ########
#########################################
function set_up_domain(L₁, L₂, N₁, N₂)

  x = Array(range(-L₁, L₁, length = 2N₁ + 1))
  y = Array(range(-L₂, L₂, length = 2N₂ + 1))

  xmid = (x[1:end-1] + x[2:end]) / 2
  ymid = (y[1:end-1] + y[2:end]) / 2

  # XMID, YMID = np.meshgrid(xmid, ymid);
  XMID = Array([x for x in xmid, y in ymid]')
  YMID = Array([y for x in xmid, y in ymid]')
  RAD = Array(sqrt.(XMID .^ 2 + YMID .^ 2))

  return XMID, YMID, RAD
end

# returns interfaces
function domain_interfaces(L₁, L₂, N₁, N₂)

  x = Array(range(-L₁, L₁, length = 2N₁ + 1))
  y = Array(range(-L₂, L₂, length = 2N₂ + 1))

  # X, Y = np.meshgrid(x, y)
  X = [a for a in x, b in y]'
  Y = [b for a in x, b in y]'
  R = sqrt.(X .^ 2 + Y .^ 2)

  X, Y, R
end

# extended domain for convoluions
function extended_domain_3_by_3(L₁, L₂, N₁, N₂)

  x = Array(range(-L₁, L₁, length = 2N₁ + 1))
  y = Array(range(-L₂, L₂, length = 2N₂ + 1))

  xext = [x .- 2 * L₁; x[2:end-1]; x .+ 2 * L₁]
  yext = [y .- 2 * L₂; y[2:end-1]; y .+ 2 * L₂]

  XEXT, YEXT = np.meshgrid(xext, yext)
  REXT = sqrt.(XEXT .^ 2 + YEXT .^ 2)

  return XEXT, YEXT, REXT
end
#########################################
#########################################
#########################################






#########################################
#### PLOTTING AND DATA VISUALISATION ####
#########################################
function plot4in1(XMID, YMID, XanD, Mel, IriD, t)
  xmid_tmp = XMID[1, :]
  N₂_tmp = Int(size(XMID)[1] / 2)

  fig = figure()
  ax = fig[:add_subplot](2, 2, 2)
  ax[:yaxis][:set_major_formatter](f)
  # subplot(221)
  title("t=$(round(t, digits=5))")
  contourf(XMID, YMID, XanD, alpha = 1, cmap = yellow, vmin = VMIN, vmax = VMAX)
  colorbar()

  ax = fig[:add_subplot](2, 2, 1)
  ax[:yaxis][:set_major_formatter](f)
  title("t=$(round(t, digits=5))")
  contourf(XMID, YMID, Mel, alpha = 1, cmap = "gray_r", vmin = VMIN, vmax = VMAX)
  colorbar()

  ax = fig[:add_subplot](2, 2, 4)
  ax[:yaxis][:set_major_formatter](f)
  contourf(XMID, YMID, IriD, alpha = 1, cmap = "Blues")
  colorbar()

  ax = fig[:add_subplot](2, 2, 3)
  ax[:yaxis][:set_major_formatter](f)
  plot(xmid_tmp, XanD[:, N₂], c = "yellow")
  plot(xmid_tmp, Mel[:, N₂], c = "black")
  plot(xmid_tmp, IriD[:, N₂], c = "blue")

  fig[:canvas][:draw]()
end

# save current plot in the format "name"_"ctr"."format"
function SavePlot(folder, name; format = "png", counter = "")
  if counter != ""
    ctrstr = @sprintf "_%07d." counter
    tmp = string(folder, name, ctrstr, format)
    savefig(tmp)
  else
    tmp = string(folder, name, ".", format)
    savefig(tmp)
  end
end

# histograms
function visualise_histogram(XINT, YINT, density; cmap)
  pcolormesh(XINT, YINT, density, cmap = cmap)
end
#########################################
#########################################
#########################################


function Convert60x60DataInto240x240Data( ABMData )
    NewData = zeros(240, 240, 151)
    for k in 1:size(ABMData, 3)
      for i in 0:240-1
          for j in 0:240-1
              NewData[i+1, j+1, k] = ABMData[floor(Int, i/4)+1, floor(Int, j/4)+1, k];
          end
      end
    end
    return NewData
end

function Convert240x240DataInto30x30Data( PDEData )
    NewData = zeros(size(PDEData, 1), 30, 30 )
    for k in 1:size(PDEData, 1)
      for i in 0:30-1
          for j in 0:30-1
              NewData[k, i+1, j+1] = mean(PDEData[k, 8*i+1:8*(i+1), 8*j+1:8*(j+1)]);
          end
      end
    end
    return NewData
end
