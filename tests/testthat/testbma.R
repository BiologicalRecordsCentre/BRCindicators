context("Test bma")

set.seed(123)

# Create data for testing
data <- data.frame(species = rep(letters, each = 50),
                   year = rep(1:50, length(letters)),
                   index = rnorm(n = 50 * length(letters), mean = 0, sd = 1), # on the unbounded scale
                   se = runif(n = 50 * length(letters), min = 0.01, max = .1))
temp <- tempfile()

test_that("simple run", {
  
  # Run the Bayesian meta-analysis
  sink(temp)
  set.seed(123)
  bma_indicator <- bma(data,
                       model = "smooth",
                       m.scale = "logit",
                       n.iter = 100,
                       seed = 123)
  sink()
  # test a data frame is returned..
  expect_is(bma_indicator, 'data.frame') 
  # it has hte right elements ...
  expect_equal(names(bma_indicator),
               c("Year", "Index.Mprime", "lowerCI.Mprime", "upperCI.Mprime", 
                 "Index.M", "lowerCI.M", "upperCI.M"))
  # and the vaues are the same
  expect_equal(bma_indicator$Index.Mprime,
               c(100, 93.3281170803379, 90.6111878405955, 87.5200429462382, 
                 85.6960730856869, 83.9068609021675, 82.5258411195666, 82.7743854797202, 
                 84.3478599908994, 86.6056218146639, 88.0782576597769, 88.5578062394473, 
                 89.238345955564, 87.982298981628, 87.8734302792336, 87.2031967135093, 
                 85.6877866411153, 84.8218459824503, 84.142244489741, 82.9412270136467, 
                 81.9500565400498, 81.4197519001476, 81.3219269630374, 82.1393227366779, 
                 83.2026871469642, 84.0002038428555, 84.4286613588619, 82.7499600249012, 
                 82.551209963727, 81.2576932647317, 80.7966839697043, 79.4798695780263, 
                 77.8596552458764, 75.5825238736814, 75.6695538589672, 75.8975875129473, 
                 77.5445898258082, 78.6391936143979, 79.8878227140963, 80.7257430563527, 
                 83.6929130770363, 85.5823427434962, 86.9888523644877, 87.5569232178518, 
                 88.3731701780635, 88.8759602490618, 90.3988588376878, 89.6151774723484, 
                 88.8136013862888, 87.5490348910981))
  expect_equal(bma_indicator$Index.M,
               c(100, 95.1066424835543, 91.3277245037226, 88.4678784823851, 
                 86.3749875308697, 84.9289578315946, 84.0282642314814, 83.5732818349891, 
                 83.4603232794704, 83.5779762787929, 83.8203922311029, 84.1037927653868, 
                 84.3668218409477, 84.5704470720489, 84.6949304070813, 84.7380799523458, 
                 84.7153660145554, 84.6562107839785, 84.5854974146199, 84.5177530126733, 
                 84.4571762815342, 84.3988976106347, 84.3321081506402, 84.2405725547339, 
                 84.1028316496622, 83.8984447989557, 83.6163659513033, 83.2555068082075, 
                 82.8247269865824, 82.3445098974237, 81.8479674471805, 81.3803882658224, 
                 80.9967396195252, 80.7479373177327, 80.6761690900314, 80.8158186574453, 
                 81.1918701375337, 81.8126218286924, 82.6678619340919, 83.7277504340206, 
                 84.9432933084802, 86.248127439241, 87.5568964352754, 88.7643412067508, 
                 89.7566668063039, 90.4240705821149, 90.6646901541542, 90.3880341473998, 
                 89.5124407395943, 87.9667954932412))
  expect_equal(bma_indicator$upperCI.Mprime,
               c(100, 101.71603026475, 100.390760343151, 97.3989075311615, 99.027915406523, 
                 93.4105829108475, 89.0609363834422, 90.6656635858461, 94.795001614849, 
                 98.7402396272953, 97.3459553696627, 97.0118402819752, 97.2850883439693, 
                 97.4541460296808, 95.677973583939, 96.8587425808644, 95.2519829741334, 
                 92.09623820932, 94.6123864016084, 92.4170748734287, 89.5064615471656, 
                 88.1619216479188, 90.5508988467181, 93.4473540316933, 92.4616676539753, 
                 92.4717453668377, 97.2948337385315, 94.1397738563561, 92.4042892449142, 
                 91.4440562370639, 90.352535255976, 90.0310366650038, 86.5736129343628, 
                 82.9753644402092, 83.9884634058, 84.4568399158517, 84.5953989600447, 
                 88.0608604261613, 88.8383050468053, 88.8849585790184, 96.4952611083576, 
                 103.638182923942, 100.144118574511, 101.545846142705, 101.524175258706, 
                 100.23781509227, 102.688628825797, 102.793775202159, 101.654403765812, 
                 106.200952061652))
  
})

test_that("degraded data", {

  data2 <- data
  # Add NAs to one year
  data2[5,c(3:4)] <- NA
  # Remove a year of data
  data2 <- data2[1:(nrow(data2)-1),]
  
  sink(temp)
  bma_indicator <- bma(data2,
                       model = "smooth",
                       m.scale = "logit",
                       n.iter = 100,
                       seed = 123)
  sink()
  
  expect_is(bma_indicator, 'data.frame')
  
})  


test_that("model options", {
  
  set.seed(123)
  sink(temp)
  bma_indicator_smooth_det2 <- bma(data,
                       model = "smooth_det2",
                       m.scale = "logit",
                       n.iter = 100,
                       seed = 123)
  sink()
  
  # test a data frame is returned..
  expect_is(bma_indicator_smooth_det2, 'data.frame')
  # it has hte right elements ...
  expect_equal(names(bma_indicator_smooth_det2),
               c("Year", "Index.Mprime", "lowerCI.Mprime", "upperCI.Mprime",
                 "Index.M", "lowerCI.M", "upperCI.M"))
  expect_equal(bma_indicator_smooth_det2$Index.M,
               c(100, 95.93909090935, 92.7672496713464, 90.3009486772026, 88.3910564873888, 
                 86.9135156978769, 85.762962294484, 84.8485307419185, 84.0907145192179, 
                 83.419528497586, 82.7801742540728, 82.1397126378509, 81.4865697868186, 
                 80.8283760531123, 80.1769336841978, 79.5395128209959, 78.918912999682, 
                 78.3155809591785, 77.7374190089791, 77.2025558362917, 76.7390450389058, 
                 76.3828846793717, 76.173517902885, 76.1536897880308, 76.3699821598276, 
                 76.8573692438152, 77.6172261446705, 78.6144721564787, 79.7775852704348, 
                 81.0389454212986, 82.367518378779, 83.7723498196311, 85.3012594012114, 
                 87.0137904449635, 88.9705565936254, 91.2341498812796, 93.8627067415033, 
                 96.8850398226239, 100.290874956475, 104.023304761472, 107.984525783209, 
                 112.054621420101, 116.091278710721, 119.928865688976, 123.379833426503, 
                 126.238820288326, 128.289739888977, 129.315919102855, 129.111127561621, 
                 127.494892855453))
  
  set.seed(123)
  
  sink(temp)
  bma_indicatorsmooth_det_sigtheta <- bma(data,
                       model = "smooth_det_sigtheta",
                       m.scale = "logit",
                       n.iter = 100,
                       seed = 123)
  sink()
  
  # test a data frame is returned..
  expect_is(bma_indicatorsmooth_det_sigtheta, 'data.frame')
  # it has the right elements ...
  expect_equal(names(bma_indicatorsmooth_det_sigtheta),
               c("Year", "Index.Mprime", "lowerCI.Mprime", "upperCI.Mprime",
                 "Index.M", "lowerCI.M", "upperCI.M"))
  expect_equal(bma_indicatorsmooth_det_sigtheta$Index.M,
               c(100, 95.1066424835543, 91.3277245037226, 88.4678784823851, 
                 86.3749875308697, 84.9289578315946, 84.0282642314814, 83.5732818349891, 
                 83.4603232794704, 83.5779762787929, 83.8203922311029, 84.1037927653868, 
                 84.3668218409477, 84.5704470720489, 84.6949304070813, 84.7380799523458, 
                 84.7153660145554, 84.6562107839785, 84.5854974146199, 84.5177530126733, 
                 84.4571762815342, 84.3988976106347, 84.3321081506402, 84.2405725547339, 
                 84.1028316496622, 83.8984447989557, 83.6163659513033, 83.2555068082075, 
                 82.8247269865824, 82.3445098974237, 81.8479674471805, 81.3803882658224, 
                 80.9967396195252, 80.7479373177327, 80.6761690900314, 80.8158186574453, 
                 81.1918701375337, 81.8126218286924, 82.6678619340919, 83.7277504340206, 
                 84.9432933084802, 86.248127439241, 87.5568964352754, 88.7643412067508, 
                 89.7566668063039, 90.4240705821149, 90.6646901541542, 90.3880341473998, 
                 89.5124407395943, 87.9667954932412))
  
})

test_that("different parameters", {
  
  sink(temp)
  set.seed(123)
  # test reading se and indexing
  bma_indicator_params1 <- bma(data,
                               model = "smooth",
                               m.scale = "logit",
                               n.iter = 100,
                               seed = 123, 
                               seFromData = TRUE,
                               rescale_indices = 2,
                               rescaleYr = 5,
                               baseline = 10,
                               incl.2deriv = TRUE)

  sink()
  
  expect_equal(bma_indicator_params1$Index.Mprime[5], 10)
  expect_equal(bma_indicator_params1$lowerCI.Mprime[5], 10)
  expect_equal(bma_indicator_params1$upperCI.Mprime[5], 10)
  
  sink(temp)
  set.seed(123)
  # test reading se
  bma_indicator_params2 <- bma(data,
                               model = "smooth",
                               m.scale = "logit",
                               n.iter = 100,
                               seed = 123, 
                               Y1perfect = FALSE)
  sink()
  expect_is(bma_indicator_params2, 'data.frame')
  
})