# A
function OuterAnnulus(r)
  return  (r .> 0.21) .* (r .< 0.25)
end

# B
function InnerBall(r)
  return (r .< 0.08)
end

# Birth Melanophores:
# need to talk about the right scaling (20 new cells per diem)
#
#     gamma_1 A convolved (XanD + IriD - beta * Mel) * (1 - Mel - XanD - IriD),
#
#where beta = 3.5 (possibly also subject to scaling)

# Death Melanophores:
# need to talk about the right scaling (#dying cells per diem)
#
#     gamma_2 A convolved (M - XanD * 1.25) * (1 - M - XanD - IriD),
#
#where beta = 3.5 (possibly also subject to scaling)



# Cell birth for IriD
#
#     gamma_3 B convolved (1 - IriD - IriL)
#
# Cell birth for IriL
#
#     gamma_3 B convolved (1 - IriD - IriL)
#
# Cell birth for XanD
#
#     gamma_3 B convolved (IriD + XanD - IriL - Mel)
#



# Switching of compartments IriD --> IriL
# think about the thresholds f,g,h not in terms of cell numbers but cell densities
#
#     A convolved (XanD / g - 1)  *  B convolved (1 - XanD/h) +  B convolved (Mel/f - 1)


# Switching of compartments IriL --> IriD
# think about the thresholds f,g,h not in terms of cell numbers but cell densities
#
#     A convolved (1 - XanD / g) * B (XanD/h - 1) + B convolved (1 - M/f)  *  B convolved (XanD/ h - 1)
