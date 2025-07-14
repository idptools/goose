# NOTE: THIS DOES NOT HAVE FRONTEND FUNCTIONALITY YET! 
# I haven't decided how I want to reimplement this yet.
# It does need to be reimplemented, however, as the previous implementation was not very user friendly.
'''
Code for desigining sequences based on epsilon. 

Epsilon is a value that is a summation of pairwise interaction values between
amino acids. It is used to predict the interaction between two proteins.
This allowes us to make sequences with epsilon values similar to a starting
sequence as well as epsilon values that are simply predicted to be overall
attractive or repulsive.  

Because epsilon is composed of the summation of an attractive and repulsive
interaction vector between two sequences, we can make sequences that have
similar interaction vectors. This can be useful for making sequences that
interact with a specific protein in a similar way to another protein. 
'''

'''
TO DO: (if time)
    • Also need seq by chunked epsilon.
'''

from collections import Counter
import random

import matplotlib.pyplot as plt
import numpy as np

from finches import epsilon_calculation
from finches.forcefields.mpipi import Mpipi_model
from finches.forcefields.calvados import calvados_model
from goose.backend import lists
from goose import goose_exceptions

# list for precomputed epsilon values
# precomputed epsilon params
precomputed_epsilon={'Mpipi_GGv1':{'A': {'A': 0.1729, 'C': 0.1479, 'D': 0.1122, 'E': 0.0991, 'F': -0.1814, 'G': 0.1355, 'H': -0.1472, 'I': 0.1184, 'K': 0.1819, 'L': 0.1302, 'M': 0.1617, 'N': 0.0295, 'P': 0.0992, 'Q': 0.0159, 'R': -0.1, 'S': 0.1578, 'T': 0.1822, 'V': 0.1273, 'W': -0.3609, 'Y': -0.2128}, 'C': {'A': 0.1479, 'C': 0.1244, 'D': 0.0871, 'E': 0.0731, 'F': -0.2188, 'G': 0.1127, 'H': -0.1831, 'I': 0.1859, 'K': 0.1586, 'L': 0.1751, 'M': 0.1461, 'N': 0.0009, 'P': 0.0673, 'Q': -0.0136, 'R': -0.1343, 'S': 0.1352, 'T': 0.16, 'V': 0.1827, 'W': -0.4051, 'Y': -0.2514}, 'D': {'A': 0.1122, 'C': 0.0871, 'D': 2.574, 'E': 2.4772, 'F': 0.0495, 'G': 0.0783, 'H': -1.1931, 'I': 0.1469, 'K': -2.5449, 'L': 0.136, 'M': 0.1068, 'N': -0.0382, 'P': 0.0172, 'Q': -0.0538, 'R': -2.4973, 'S': 0.0989, 'T': 0.1225, 'V': 0.1447, 'W': -0.0261, 'Y': 0.0362}, 'E': {'A': 0.0991, 'C': 0.0731, 'D': 2.4772, 'E': 2.3835, 'F': 0.0398, 'G': 0.0643, 'H': -1.1415, 'I': 0.1344, 'K': -2.4518, 'L': 0.1232, 'M': 0.0931, 'N': -0.056, 'P': -0.0018, 'Q': -0.0721, 'R': -2.3995, 'S': 0.0854, 'T': 0.1096, 'V': 0.1323, 'W': -0.0378, 'Y': 0.0262}, 'F': {'A': -0.1814, 'C': -0.2188, 'D': 0.0495, 'E': 0.0398, 'F': -0.6114, 'G': -0.2051, 'H': -0.566, 'I': -0.1716, 'K': 0.0014, 'L': -0.1832, 'M': -0.2144, 'N': -0.3574, 'P': -0.393, 'Q': -0.3819, 'R': -0.2339, 'S': -0.199, 'T': -0.1842, 'V': -0.1662, 'W': -0.8217, 'Y': -0.6489}, 'G': {'A': 0.1355, 'C': 0.1127, 'D': 0.0783, 'E': 0.0643, 'F': -0.2051, 'G': -0.0747, 'H': -0.1708, 'I': 0.1667, 'K': 0.1411, 'L': 0.1568, 'M': 0.1302, 'N': -0.0007, 'P': 0.0525, 'Q': -0.0154, 'R': -0.1286, 'S': -0.0354, 'T': 0.1446, 'V': 0.1647, 'W': -0.3788, 'Y': -0.2356}, 'H': {'A': -0.1472, 'C': -0.1831, 'D': -1.1931, 'E': -1.1415, 'F': -0.566, 'G': -0.1708, 'H': 0.8857, 'I': -0.1362, 'K': 1.3994, 'L': -0.1475, 'M': -0.178, 'N': -0.3183, 'P': -0.3451, 'Q': -0.3419, 'R': 1.289, 'S': -0.1641, 'T': -0.1492, 'V': -0.1313, 'W': -0.7715, 'Y': -0.6027}, 'I': {'A': 0.1184, 'C': 0.1859, 'D': 0.1469, 'E': 0.1344, 'F': -0.1716, 'G': 0.1667, 'H': -0.1362, 'I': 0.0416, 'K': 0.2281, 'L': 0.0644, 'M': 0.1023, 'N': 0.056, 'P': 0.149, 'Q': 0.0432, 'R': -0.0806, 'S': 0.1954, 'T': 0.2247, 'V': 0.0536, 'W': -0.3651, 'Y': -0.2054}, 'K': {'A': 0.1819, 'C': 0.1586, 'D': -2.5449, 'E': -2.4518, 'F': 0.0014, 'G': 0.1411, 'H': 1.3994, 'I': 0.2281, 'K': 2.2789, 'L': 0.2165, 'M': 0.1853, 'N': 0.0263, 'P': 0.1123, 'Q': 0.0126, 'R': 2.0644, 'S': 0.1688, 'T': 0.1977, 'V': 0.2233, 'W': 0.017, 'Y': 0.0275}, 'L': {'A': 0.1302, 'C': 0.1751, 'D': 0.136, 'E': 0.1232, 'F': -0.1832, 'G': 0.1568, 'H': -0.1475, 'I': 0.0644, 'K': 0.2165, 'L': 0.0794, 'M': 0.1153, 'N': 0.045, 'P': 0.1346, 'Q': 0.0319, 'R': -0.0924, 'S': 0.1849, 'T': 0.2138, 'V': 0.0761, 'W': -0.3771, 'Y': -0.2171}, 'M': {'A': 0.1617, 'C': 0.1461, 'D': 0.1068, 'E': 0.0931, 'F': -0.2144, 'G': 0.1302, 'H': -0.178, 'I': 0.1023, 'K': 0.1853, 'L': 0.1153, 'M': 0.1501, 'N': 0.0155, 'P': 0.0957, 'Q': 0.0016, 'R': -0.1241, 'S': 0.1566, 'T': 0.1844, 'V': 0.1118, 'W': -0.4093, 'Y': -0.2485}, 'N': {'A': 0.0295, 'C': 0.0009, 'D': -0.0382, 'E': -0.056, 'F': -0.3574, 'G': -0.0007, 'H': -0.3183, 'I': 0.056, 'K': 0.0263, 'L': 0.045, 'M': 0.0155, 'N': -0.127, 'P': -0.0983, 'Q': -0.1452, 'R': -0.2736, 'S': 0.0151, 'T': 0.0356, 'V': 0.0562, 'W': -0.5511, 'Y': -0.3917}, 'P': {'A': 0.0992, 'C': 0.0673, 'D': 0.0172, 'E': -0.0018, 'F': -0.393, 'G': 0.0525, 'H': -0.3451, 'I': 0.149, 'K': 0.1123, 'L': 0.1346, 'M': 0.0957, 'N': -0.0983, 'P': 0.0539, 'Q': -0.1179, 'R': -0.2801, 'S': 0.0821, 'T': 0.1148, 'V': 0.1451, 'W': -0.6428, 'Y': -0.4368}, 'Q': {'A': 0.0159, 'C': -0.0136, 'D': -0.0538, 'E': -0.0721, 'F': -0.3819, 'G': -0.0154, 'H': -0.3419, 'I': 0.0432, 'K': 0.0126, 'L': 0.0319, 'M': 0.0016, 'N': -0.1452, 'P': -0.1179, 'Q': -0.1638, 'R': -0.2956, 'S': 0.001, 'T': 0.0222, 'V': 0.0435, 'W': -0.5805, 'Y': -0.417}, 'R': {'A': -0.1, 'C': -0.1343, 'D': -2.4973, 'E': -2.3995, 'F': -0.2339, 'G': -0.1286, 'H': 1.289, 'I': -0.0806, 'K': 2.0644, 'L': -0.0924, 'M': -0.1241, 'N': -0.2736, 'P': -0.2801, 'Q': -0.2956, 'R': 2.08, 'S': -0.1167, 'T': -0.098, 'V': -0.0777, 'W': -0.3863, 'Y': -0.2948}, 'S': {'A': 0.1578, 'C': 0.1352, 'D': 0.0989, 'E': 0.0854, 'F': -0.199, 'G': -0.0354, 'H': -0.1641, 'I': 0.1954, 'K': 0.1688, 'L': 0.1849, 'M': 0.1566, 'N': 0.0151, 'P': 0.0821, 'Q': 0.001, 'R': -0.1167, 'S': 0.1455, 'T': 0.1698, 'V': 0.1921, 'W': -0.3809, 'Y': -0.2308}, 'T': {'A': 0.1822, 'C': 0.16, 'D': 0.1225, 'E': 0.1096, 'F': -0.1842, 'G': 0.1446, 'H': -0.1492, 'I': 0.2247, 'K': 0.1977, 'L': 0.2138, 'M': 0.1844, 'N': 0.0356, 'P': 0.1148, 'Q': 0.0222, 'R': -0.098, 'S': 0.1698, 'T': 0.1965, 'V': 0.2204, 'W': -0.371, 'Y': -0.2168}, 'V': {'A': 0.1273, 'C': 0.1827, 'D': 0.1447, 'E': 0.1323, 'F': -0.1662, 'G': 0.1647, 'H': -0.1313, 'I': 0.0536, 'K': 0.2233, 'L': 0.0761, 'M': 0.1118, 'N': 0.0562, 'P': 0.1451, 'Q': 0.0435, 'R': -0.0777, 'S': 0.1921, 'T': 0.2204, 'V': 0.0726, 'W': -0.3555, 'Y': -0.1992}, 'W': {'A': -0.3609, 'C': -0.4051, 'D': -0.0261, 'E': -0.0378, 'F': -0.8217, 'G': -0.3788, 'H': -0.7715, 'I': -0.3651, 'K': 0.017, 'L': -0.3771, 'M': -0.4093, 'N': -0.5511, 'P': -0.6428, 'Q': -0.5805, 'R': -0.3863, 'S': -0.3809, 'T': -0.371, 'V': -0.3555, 'W': -1.0436, 'Y': -0.8617}, 'Y': {'A': -0.2128, 'C': -0.2514, 'D': 0.0362, 'E': 0.0262, 'F': -0.6489, 'G': -0.2356, 'H': -0.6027, 'I': -0.2054, 'K': 0.0275, 'L': -0.2171, 'M': -0.2485, 'N': -0.3917, 'P': -0.4368, 'Q': -0.417, 'R': -0.2948, 'S': -0.2308, 'T': -0.2168, 'V': -0.1992, 'W': -0.8617, 'Y': -0.687}},
                     'CALVADOS2':{'A': {'A': 0.5025, 'C': 0.1929, 'D': 0.7216, 'E': 0.7598, 'F': -0.1997, 'G': 0.123, 'H': 0.257, 'I': 0.1691, 'K': 0.5583, 'L': 0.0588, 'M': 0.1816, 'N': 0.3213, 'P': 0.3955, 'Q': 0.3383, 'R': -0.0649, 'S': 0.3105, 'T': 0.3804, 'V': 0.5413, 'W': -0.3764, 'Y': -0.3311}, 'C': {'A': 0.1929, 'C': -0.14, 'D': 0.4074, 'E': 0.4384, 'F': -0.5668, 'G': -0.1903, 'H': -0.0873, 'I': -0.1807, 'K': 0.2186, 'L': -0.295, 'M': -0.1677, 'N': -0.0112, 'P': 0.0689, 'Q': -0.0015, 'R': -0.4314, 'S': -0.0105, 'T': 0.0517, 'V': 0.213, 'W': -0.7582, 'Y': -0.705}, 'D': {'A': 0.7216, 'C': 0.4074, 'D': 2.572, 'E': 2.5294, 'F': 0.0198, 'G': 0.3079, 'H': 0.2894, 'I': 0.4004, 'K': -4.2734, 'L': 0.2851, 'M': 0.4134, 'N': 0.5471, 'P': 0.622, 'Q': 0.5733, 'R': -4.816, 'S': 0.5233, 'T': 0.6076, 'V': 0.7819, 'W': -0.1531, 'Y': -0.1147}, 'E': {'A': 0.7598, 'C': 0.4384, 'D': 2.5294, 'E': 2.4898, 'F': 0.0454, 'G': 0.33, 'H': 0.3328, 'I': 0.435, 'K': -4.0533, 'L': 0.3166, 'M': 0.4484, 'N': 0.5831, 'P': 0.6596, 'Q': 0.6118, 'R': -4.6153, 'S': 0.556, 'T': 0.6451, 'V': 0.8254, 'W': -0.129, 'Y': -0.0918}, 'F': {'A': -0.1997, 'C': -0.5668, 'D': 0.0198, 'E': 0.0454, 'F': -1.0373, 'G': -0.6043, 'H': -0.5206, 'I': -0.6222, 'K': -0.1988, 'L': -0.7444, 'M': -0.6083, 'N': -0.432, 'P': -0.3437, 'Q': -0.4277, 'R': -0.8954, 'S': -0.422, 'T': -0.3633, 'V': -0.1948, 'W': -1.2459, 'Y': -1.1859}, 'G': {'A': 0.123, 'C': -0.1903, 'D': 0.3079, 'E': 0.33, 'F': -0.6043, 'G': -0.2117, 'H': -0.156, 'I': -0.2442, 'K': 0.1195, 'L': -0.3494, 'M': -0.2323, 'N': -0.0771, 'P': -0.0013, 'Q': -0.0758, 'R': -0.4846, 'S': -0.0652, 'T': -0.0183, 'V': 0.1245, 'W': -0.7926, 'Y': -0.7344}, 'H': {'A': 0.257, 'C': -0.0873, 'D': 0.2894, 'E': 0.3328, 'F': -0.5206, 'G': -0.156, 'H': -0.0156, 'I': -0.1193, 'K': 0.4863, 'L': -0.239, 'M': -0.1057, 'N': 0.0511, 'P': 0.1335, 'Q': 0.0662, 'R': -0.1943, 'S': 0.0445, 'T': 0.1164, 'V': 0.289, 'W': -0.7131, 'Y': -0.6634}, 'I': {'A': 0.1691, 'C': -0.1807, 'D': 0.4004, 'E': 0.435, 'F': -0.6222, 'G': -0.2442, 'H': -0.1193, 'I': -0.2169, 'K': 0.2056, 'L': -0.3376, 'M': -0.2033, 'N': -0.0424, 'P': 0.0414, 'Q': -0.0292, 'R': -0.4771, 'S': -0.0461, 'T': 0.0237, 'V': 0.1963, 'W': -0.8181, 'Y': -0.7665}, 'K': {'A': 0.5583, 'C': 0.2186, 'D': -4.2734, 'E': -4.0533, 'F': -0.1988, 'G': 0.1195, 'H': 0.4863, 'I': 0.2056, 'K': 2.2306, 'L': 0.0834, 'M': 0.2195, 'N': 0.3657, 'P': 0.4466, 'Q': 0.3905, 'R': 1.9968, 'S': 0.3448, 'T': 0.4307, 'V': 0.6138, 'W': -0.3831, 'Y': -0.3416}, 'L': {'A': 0.0588, 'C': -0.295, 'D': 0.2851, 'E': 0.3166, 'F': -0.7444, 'G': -0.3494, 'H': -0.239, 'I': -0.3376, 'K': 0.0834, 'L': -0.4583, 'M': -0.3239, 'N': -0.1585, 'P': -0.0736, 'Q': -0.1484, 'R': -0.6011, 'S': -0.1576, 'T': -0.0919, 'V': 0.0785, 'W': -0.944, 'Y': -0.8896}, 'M': {'A': 0.1816, 'C': -0.1677, 'D': 0.4134, 'E': 0.4484, 'F': -0.6083, 'G': -0.2323, 'H': -0.1057, 'I': -0.2033, 'K': 0.2195, 'L': -0.3239, 'M': -0.1896, 'N': -0.0292, 'P': 0.0545, 'Q': -0.0157, 'R': -0.4631, 'S': -0.0334, 'T': 0.0368, 'V': 0.2097, 'W': -0.8039, 'Y': -0.7525}, 'N': {'A': 0.3213, 'C': -0.0112, 'D': 0.5471, 'E': 0.5831, 'F': -0.432, 'G': -0.0771, 'H': 0.0511, 'I': -0.0424, 'K': 0.3657, 'L': -0.1585, 'M': -0.0292, 'N': 0.1227, 'P': 0.2024, 'Q': 0.1374, 'R': -0.2915, 'S': 0.1161, 'T': 0.1858, 'V': 0.3532, 'W': -0.6199, 'Y': -0.5708}, 'P': {'A': 0.3955, 'C': 0.0689, 'D': 0.622, 'E': 0.6596, 'F': -0.3437, 'G': -0.0013, 'H': 0.1335, 'I': 0.0414, 'K': 0.4466, 'L': -0.0736, 'M': 0.0545, 'N': 0.2024, 'P': 0.2807, 'Q': 0.2187, 'R': -0.2036, 'S': 0.1933, 'T': 0.2646, 'V': 0.4316, 'W': -0.5282, 'Y': -0.4808}, 'Q': {'A': 0.3383, 'C': -0.0015, 'D': 0.5733, 'E': 0.6118, 'F': -0.4277, 'G': -0.0758, 'H': 0.0662, 'I': -0.0292, 'K': 0.3905, 'L': -0.1484, 'M': -0.0157, 'N': 0.1374, 'P': 0.2187, 'Q': 0.1545, 'R': -0.2822, 'S': 0.1278, 'T': 0.202, 'V': 0.3753, 'W': -0.6172, 'Y': -0.5693}, 'R': {'A': -0.0649, 'C': -0.4314, 'D': -4.816, 'E': -4.6153, 'F': -0.8954, 'G': -0.4846, 'H': -0.1943, 'I': -0.4771, 'K': 1.9968, 'L': -0.6011, 'M': -0.4631, 'N': -0.2915, 'P': -0.2036, 'Q': -0.2822, 'R': 1.7581, 'S': -0.2888, 'T': -0.2226, 'V': -0.0482, 'W': -1.1003, 'Y': -1.0445}, 'S': {'A': 0.3105, 'C': -0.0105, 'D': 0.5233, 'E': 0.556, 'F': -0.422, 'G': -0.0652, 'H': 0.0445, 'I': -0.0461, 'K': 0.3448, 'L': -0.1576, 'M': -0.0334, 'N': 0.1161, 'P': 0.1933, 'Q': 0.1278, 'R': -0.2888, 'S': 0.1136, 'T': 0.177, 'V': 0.3359, 'W': -0.6071, 'Y': -0.5565}, 'T': {'A': 0.3804, 'C': 0.0517, 'D': 0.6076, 'E': 0.6451, 'F': -0.3633, 'G': -0.0183, 'H': 0.1164, 'I': 0.0237, 'K': 0.4307, 'L': -0.0919, 'M': 0.0368, 'N': 0.1858, 'P': 0.2646, 'Q': 0.202, 'R': -0.2226, 'S': 0.177, 'T': 0.2484, 'V': 0.416, 'W': -0.5488, 'Y': -0.5011}, 'V': {'A': 0.5413, 'C': 0.213, 'D': 0.7819, 'E': 0.8254, 'F': -0.1948, 'G': 0.1245, 'H': 0.289, 'I': 0.1963, 'K': 0.6138, 'L': 0.0785, 'M': 0.2097, 'N': 0.3532, 'P': 0.4316, 'Q': 0.3753, 'R': -0.0482, 'S': 0.3359, 'T': 0.416, 'V': 0.5909, 'W': -0.3763, 'Y': -0.3334}, 'W': {'A': -0.3764, 'C': -0.7582, 'D': -0.1531, 'E': -0.129, 'F': -1.2459, 'G': -0.7926, 'H': -0.7131, 'I': -0.8181, 'K': -0.3831, 'L': -0.944, 'M': -0.8039, 'N': -0.6199, 'P': -0.5282, 'Q': -0.6172, 'R': -1.1003, 'S': -0.6071, 'T': -0.5488, 'V': -0.3763, 'W': -1.4609, 'Y': -1.399}, 'Y': {'A': -0.3311, 'C': -0.705, 'D': -0.1147, 'E': -0.0918, 'F': -1.1859, 'G': -0.7344, 'H': -0.6634, 'I': -0.7665, 'K': -0.3416, 'L': -0.8896, 'M': -0.7525, 'N': -0.5708, 'P': -0.4808, 'Q': -0.5693, 'R': -1.0445, 'S': -0.5565, 'T': -0.5011, 'V': -0.3334, 'W': -1.399, 'Y': -1.3363}}}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#        -=-=-=-=-=-= Code for epsilon related stuffs =-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def load_IMC_object(model='Mpipi_GGv1'):
    '''
    Function to load a FINCHES IMC object.
    Lets us only load it once and then can use it iteratively. 

    Parameters
    -----------
    model : string
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    '''

    # check the current implementations of forcefields.
    if model not in ['Mpipi_GGv1', 'CALVADOS2']:
        raise goose_exceptions.GooseInputError('Only Mpipi_GGv1 and CALVADOS2 forcefields have been implemented.')
    
    # initialize forcefield parameters
    if model in ['Mpipi_GGv1']:
        ff_model = Mpipi_model(model)
    else:
        ff_model = calvados_model(model)
    
    # make IMC Object
    IMC_object = epsilon_calculation.InteractionMatrixConstructor(ff_model)
    return IMC_object

def get_interaction_vectors(seq1, seq2, IMC_object, approach='mean'):
    '''
    Function that returns the interaction vectors between two sequences. 
    The vector will return a list where the first list is seq1 vector
    and the second item in the list is seq2 vector. 

    Should be able to use this to map reigions in our protein that need modifications
    to comply with an objective interaction vector. 

    Parameters
    -----------
    seq1 : string
        the amino acid sequence for seq1 as a string

    seq2 : string
        the amino acid sequence for seq2 as a string

    IMC_object : Object
        A loaded FINCHES IMC object.

    approach : string
        How the vectors should be calculated. Options are sum and mean.
        Default = sum. 
        Mean (should) give approximate values per amino acid in the sequence which
        in theory will identify interaction propensities that are closer to values
        in the pairwise interaction dict. 

    Returns
    -------
    list of np.arrays
        Returns a list of np.arrays. First item in the list is the
        interaction vector for seq1, the second is for seq2. 
    '''
    # use IMC_object to calculate heterotypic matrix
    interaction_matrix=IMC_object.calculate_pairwise_heterotypic_matrix(seq1,seq2)
    # now get the mean or sum values
    if approach=='mean':
        t1_vector=(np.mean(interaction_matrix, axis=1))
        t2_vector=(np.mean(interaction_matrix, axis=0))
    elif approach == 'sum':
        t1_vector=(np.sum(interaction_matrix, axis=1))
        t2_vector=(np.sum(interaction_matrix, axis=0))
    else:
        raise Exception('approach can only be set to mean or sum at the moment.')

    return [t1_vector, t2_vector]


def get_epsilon_vectors(seq1, seq2, IMC_object):
    '''
    Function that returns the interaction vectors between two sequences. 
    The vector will return a list where the first list is seq1 vector
    and the second item in the list is seq2 vector. 

    Should be able to use this to map reigions in our protein that need modifications
    to comply with an objective interaction vector. 

    Parameters
    -----------
    seq1 : string
        the amino acid sequence for seq1 as a string

    seq2 : string
        the amino acid sequence for seq2 as a string

    IMC_object : Object
        A loaded FINCHES IMC object.

    Returns
    -------
    list of np.arrays
        Returns a list of np.arrays. First item in the list is the
        interaction vector for seq1, the second is for seq2. 
    '''
    # use IMC_object to calculate heterotypic matrix
    interaction_matrix=IMC_object.calculate_epsilon_vectors(seq1,seq2)

    return [interaction_matrix[0], interaction_matrix[1]]


def get_epsilon_value(seq1, seq2, IMC_object):
    '''
    Function to get the epsilon value using the Mpipi_GGv1 model. 

    Parameters
    -----------
    seq1 : string
        the amino acid sequence for seq1 as a string

    seq2 : string
        the amino acid sequence for seq2 as a string

    IMC_object : Object
        A loaded FINCHES IMC object.

    Returns
    -------
    list of np.arrays
        Returns a list of np.arrays. First item in the list is the
        interaction vector for seq1, the second is for seq2. 
    '''    
    # use IMC_object to calculate heterotypic matrix
    epslon_value=IMC_object.calculate_epsilon_value(seq1,seq2)

    return epslon_value


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#    -=-=-=-=-=-=-=-=- Code for modifiying sequences =-=-=-=-=-=-=-=-=-
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



def optimize_to_epsilon_value(starting_sequence, interacting_sequence, objective_epsilon,
    allowed_error=0.1, optimization_iterations=None, exclude_aas=[], model='Mpipi_GGv1',
    preloaded_IMC_object=None, return_best_sequence=False, maximal_optimization=False):
    '''
    Function to optimze a sequence to have a specific epsilon value relative to another sequence. 

    Parameters
    -----------
    starting_sequence : str
        The sequence that you want to optimize to have a specific epsilon value relative to
        the interacting sequence.

    interacting_sequence : str
        The sequence that the starting sequence is interacting with.

    objective_epsilon : float
        The epsilon value that you want the starting sequence to have relative to the interacting

    allowed_error : float
        The allowed error in the epsilon value. Default is 0.1

    optimization_iterations : int
        The number of iterations to run the optimization. Default is length of the sequence

    exclude_aas : list
        A list of amino acids to exclude from being added to the sequence during optimization.

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    preloaded_IMC_object : object
        If we are using this optimizer and we already loaded the FINCHES model, we can
        skip that here to speed things up. 
        Default is None. 

    return_best_sequence : bool
        whether to just return the best sequence we could get to
        the objective epsilon value. 
        Default is False. 

    maximal_optimization : bool
        Whether to optimize to the maximal extent possible. 
        Reduces the sequence space explored, so default is False

    Returns
    -------
    str
        The optimized sequence that has the closest epsilon value to the objective epsilon value.
    '''

    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas=precomputed_epsilon[model]
    
    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
    # make list of exclude aas non redundant
    exclude_aas=list(set(exclude_aas))

    # make sure we aren't removing all of the amino acids
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all amino acids.')

    # get rid of amino acids as needed
    for aa in exclude_aas:
        all_aas.remove(aa)

    # if optimization_iterations=-None, set to lenght of sequence
    if optimization_iterations is None:
        optimization_iterations=len(starting_sequence)

    if preloaded_IMC_object==None:
        # load IMC object
        loaded_model = load_IMC_object(model)
    else:
        loaded_model=preloaded_IMC_object

    # calculate the epsilon value of the starting sequence
    eps_vectors=get_epsilon_vectors(starting_sequence, interacting_sequence, loaded_model)
    starting_epsilon=np.sum(eps_vectors)
    # calculate the difference between the starting epsilon and the objective epsilon
    epsilon_diff=objective_epsilon-starting_epsilon
    # set the epsilon difference to be the current difference
    current_epsilon_diff=epsilon_diff
    # set the current sequence to be the starting sequence
    current_sequence=starting_sequence

    # now iterate over sequence to get towards the objective value. 
    for i in range(0, optimization_iterations):
        # if i > 1, recalculate the eps vectors. 
        if i > 1:
            eps_vectors=get_epsilon_vectors(current_sequence, interacting_sequence, loaded_model)
        # choose a random amino acid to change. 
        aa_to_change=random.randint(0, len(starting_sequence)-1)
        
        # dict to hold possible aas
        possible_aas={}

        # iterate through all aas
        for aa in all_aas:
            curv=0
            for a in interacting_sequence:
                curv=curv+epsilon_aas[aa][a]
            possible_aas[aa]=curv/len(interacting_sequence)
        # sort the possible aas by the difference between the current epsilon value and the objective epsilon value
        possible_aas=sorted(possible_aas, key=lambda x:abs(possible_aas[x]-current_epsilon_diff))

        # see if we are trying to maximize the optimziation
        if maximal_optimization==True:
            new_aa=possible_aas[0]
        else:
            # based on the class, allow variable amounts of 'randomness'
            # this lets us keep the chemical interaction specificity close
            # but explore maximal sequecne space.
            if possible_aas[0] in ['D', 'E', 'K', 'R']:
                randomness=2
            elif possible_aas[0] in ['W', 'Y', 'F']:
                randomness=3
            elif possible_aas[0] in ['A', 'I', 'L', 'V', 'M']:
                randomness=5
            elif possible_aas[0] in ['N', 'Q', 'S', 'T', 'G']:
                randomness=4
            else:
                randomness=3
            # choose amino acid to add
            new_aa=random.choice(possible_aas[:randomness])


        # change the amino acid
        current_sequence=current_sequence[:aa_to_change]+new_aa+current_sequence[aa_to_change+1:]
        # recalculate the epsilon value
        eps_vectors=get_epsilon_vectors(current_sequence, interacting_sequence, loaded_model)
        current_epsilon=np.sum(eps_vectors)
        # calculate the new difference
        current_epsilon_diff=objective_epsilon-current_epsilon
        # if the difference is less than the allowed error, break
        if abs(current_epsilon_diff) < allowed_error:
            return current_sequence
    if return_best_sequence==False:
        raise goose_exceptions.GooseFail('Failed to optimize to epsilon value.')
    else:
        return current_sequence



def return_constrained_aa_list(input_sequence, exclude_aas=[]):
    '''
    returns a constrained aa list depending on if we get too many of
    any specific amino acid. 
    
    Parameters
    -----------
    input_sequence : str
        The sequence that we are checking for amino acid counts

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence.

    Returns
    -------
    list
        A list of amino acids that are allowed to be used in the sequence

    '''
    # list all aas
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    if input_sequence=='':
        finaas=[]
        for aa in all_aas:
            if aa not in exclude_aas:
                finaas.append(aa)
        return finaas

    # lists we will watch out for 
    aro=['W', 'F', 'Y']
    ali=['I', 'L', 'V', 'A','M']

    # make list of exclude aas non redundant
    exclude_aas=list(set(exclude_aas))

    # make sure not all aas are excluded
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all amino acids.')

    # exclude some aa if needed
    for aa in exclude_aas:
        all_aas.remove(aa)

    # get the count of each amino acid
    aa_counts=Counter(input_sequence)
    
    total_aro = sum([aa_counts[aa] for aa in aro])
    total_ali = sum([aa_counts[aa] for aa in ali])
    
    # make sure we don't get a weird amount of C
    if 'C' in all_aas:
        if 'C' in input_sequence:
            if aa_counts['C']/len(input_sequence) > 0.15:
                all_aas.remove('C')

    # make sure we don't get too many aro
    if total_aro/len(input_sequence) > 0.2:
        for aa in aro:
            # check in case we already excluded it
            if aa in all_aas:
                all_aas.remove(aa)

    # modulate maximum allowable aliphatics dependent on aromaitcs. 
    max_ali = 0.3 - total_aro/len(input_sequence)

    # make sure we don't get too many ali
    if total_ali/len(input_sequence) > max_ali:
        for aa in ali:
            if aa in all_aas:
                all_aas.remove(aa)

    # make sure all aas not equal to []
    if all_aas == []:
        # if list is empty, we will add anything the user didn't explicitly exclude. 
        # this is very unlikely, but we should priortize that over the fractions.
        for aa in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
            if aa not in exclude_aas:
                all_aas.append(aa)

    return all_aas

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#    -=-=-=-=-=-=-=-=-=-=- Code for sequence creation =-=-=-=-=-=-=-=-=-
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def create_starting_sequence(length, exclude_aas=[]):
    '''
    Function to create a random amino acid sequence of a specific length.
    The randomness is contstrained at least to increase probability
    of making something disordered but otherwise generally doesn't make
    the same sequence twice.

    Parameters
    -----------
    length : int
        The length of the sequence you want to create

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence.
        Default = []

    Returns
    -------
    str
        The random amino acid sequence of the specified length
    '''
    if length <= 0:
        raise goose_exceptions.GooseInputError('Length must be greater than 0.')

    # make list of exclude aas non redundant
    exclude_aas=list(set(exclude_aas))

    if exclude_aas != []:
        if len(exclude_aas)==20:
            raise goose_exceptions.GooseInputError('Cannot exclude all amino acids.')
        fl=[]
        for aa in lists.disordered_list:
            if aa not in exclude_aas:
                fl.append(aa)
    else:
        fl=lists.disordered_list

    return ''.join(np.random.choice(fl, length))


def create_seq_by_epsilon_vectors(sequence_of_interest, interacting_sequence=None, 
    exclude_aas=[], model='Mpipi_GGv1'):
    '''
    Function to create a sequence with approximately similar 
    interaction vectors to a sequence of interest. This function
    has some tolerance in that it will randomly choose from the
    best 2 amino acids at each position, so the total epsilon value
    will differ. However, the patterna cross the sequence should be quite close. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want your generated sequence to have an epsilon value relative to.

    interacting_sequence : str
        If you want to make a sequence that has an epsilon equal
        to that between your sequence_of_interest and some interacting sequence,
        set this equal to the interacting sequence. If None, then it will be the same as sequence of interest (
        so will make something with the same homotypic interactions).

        ex. If you want to make a FUS variant that interacts with hnRNPA1 in a specific way, 
        set sequence_of_interest=FUS and interacting_sequence=hnRNPA1.

        Default = None

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas=precomputed_epsilon[model]
    
    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all amino acids.')

    # remove AAs we are excluding
    for a in exclude_aas:
        all_aas.remove(a)

    # see if we need to make interacting sequence the same as sequence of interest
    if interacting_sequence is None:
        interacting_sequence = sequence_of_interest

    # load imc object
    loaded_model = load_IMC_object(model)

    # calculate eps vectors
    eps_vectors=get_epsilon_vectors(sequence_of_interest, interacting_sequence, loaded_model)
    # sum them to get overall interaction vectors
    interaction_vectors=eps_vectors[0]+eps_vectors[1]

    # make empty string to hold new sequence
    new_sequence = ''

    # now we need to choose amino acids that match the interaction vector of our 
    # sequence of interest. 
    for i in range(0, len(sequence_of_interest)):
        # get the interaction value for the amino acid
        aa_interaction=interaction_vectors[i]
        # get the amino acids that are closest to the interaction value
        possible_vals={}
        for aa in all_aas:
            cur_tot=0
            for aa_interacting in interacting_sequence:
                cur_tot=cur_tot+epsilon_aas[aa][aa_interacting]
            possible_vals[aa]=cur_tot/len(interacting_sequence)
        possible_aas=sorted(possible_vals, key=lambda x:abs(possible_vals[x]-aa_interaction))
        # based on the class, allow variable amounts of 'randomness'
        # this lets us keep the chemical interaction specificity close
        # but explore maximal sequecne space.
        if possible_aas[0] in ['D', 'E', 'K', 'R']:
            randomness=2
        elif possible_aas[0] in ['W', 'Y', 'F']:
            randomness=3
        elif possible_aas[0] in ['A', 'I', 'L', 'V', 'M']:
            randomness=5
        elif possible_aas[0] in ['N', 'Q', 'S', 'T', 'G']:
            randomness=4
        else:
            randomness=3
        # choose amino acid to add
        new_sequence=new_sequence+random.choice(possible_aas[:randomness])
    return new_sequence



def create_seq_by_epsilon(sequence_of_interest, interacting_sequence=None, 
    exclude_aas=[], model='Mpipi_GGv1'):
    '''
    Function that creates a sequence with the same total matching epsilon. Doesn't
    Care about linear space as much as previous function. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want your generated sequence to have an epsilon value relative to.

    interacting_sequence : str
        If you want to make a sequence that has an epsilon equal
        to that between your sequence_of_interest and some interacting sequence,
        set this equal to the interacting sequence. If None, then it will be the same as sequence of interest (
        so will make something with the same homotypic interactions).
        
        ex. If you want to make a FUS variant that interacts with hnRNPA1 in a specific way, 
        set sequence_of_interest=FUS and interacting_sequence=hnRNPA1.

        Default = None

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []        

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas=precomputed_epsilon[model]
    
    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all aminon acids.')

    # remove excluded amino acids
    for aa in exclude_aas:
        all_aas.remove(aa)

    # make sequence interactor self if None is specified.
    if interacting_sequence==None:
        interacting_sequence=sequence_of_interest
    
    # load imc object
    loaded_model = load_IMC_object(model)

    # get epsilon value between sequence and interacting_sequence
    eps_vectors=get_epsilon_vectors(sequence_of_interest, interacting_sequence, loaded_model)
    epsilon = np.sum(eps_vectors)
    # make a starting sequence equal to length of the starting sequence. 
    starting_sequence=create_starting_sequence(len(sequence_of_interest), exclude_aas=exclude_aas)
    # using the random sequence as a starting point, now optimize back to the epsilon value
    # we want
    new_sequence=optimize_to_epsilon_value(starting_sequence, interacting_sequence, 
        epsilon, exclude_aas=exclude_aas, model=model, preloaded_IMC_object=loaded_model)
    return new_sequence



def create_attractive_or_repulsive_seq(objective_seq_length, interacting_sequence, 
    attractive_or_repulsive, exclude_aas=[], model='Mpipi_GGv1'):
    '''
    Function that creates a sequence that is attractive
    or repulsive to the interacting sequence. You can also
    specify the length of the sequence. 
    
    Parameters
    -----------
    objective_seq_length : int
        The length of the sequence that you want to generate.

    interacting_sequence : str
        The sequence that the sequence of interest is interacting with.

    attractive_or_repulsive : str
        Whether the sequence should be attractive or repulsive to the interacting sequence.

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # load the model
    loaded_model = load_IMC_object(model)

    if attractive_or_repulsive not in ['attractive', 'repulsive']:
        raise goose_exceptions.GooseInputError("attractive_or_repulsive must be either 'attractive' or 'repulsive'.")

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas=precomputed_epsilon[model]

    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    # check number of amino acids excluded
    if len(set(exclude_aas))==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all aminon acids.')


    # string to hold new sequence we are generating
    new_sequence=''

    # loop through the length of the sequence we want to generate
    for i in range(0, objective_seq_length):
        # hold possible aas
        possible_aas={}

        # if seq over 8 amino acids, start constrainin all_aas
        if len(new_sequence) > 8:
            all_aas=return_constrained_aa_list(new_sequence, exclude_aas=exclude_aas)
        else:
            all_aas=return_constrained_aa_list('', exclude_aas=exclude_aas)

        # loop through all aas
        for aa in all_aas:
            total_interaction=0
            # loop through all aas in the interacting sequence
            for aa2 in interacting_sequence:
                # get the epsilon value between the aa and the interacting sequence
                interaction_value=total_interaction+epsilon_aas[aa][aa2]
            possible_aas[aa]=interaction_value
        # choose all amino acids with + or - value depending on if we want attractive or repulsive
        if attractive_or_repulsive=='attractive':
            use_these_aas={k:v for k,v in possible_aas.items() if v<0}
        else:
            use_these_aas={k:v for k,v in possible_aas.items() if v>0}

        # if we have fewer than 3 amino acids...
        if len(use_these_aas) < 3:
            if attractive_or_repulsive=='attractive':
                # sort by most negative to positive
                sorted_aas=sorted(possible_aas.items(), key=lambda x: x[1])
            else:
                # sort by most negative to positive
                sorted_aas=sorted(possible_aas.items(), key=lambda x: x[1], reverse=True)
            # add from best 5
            new_sequence+=random.choice(sorted_aas[:3])[0]
        else:
            if attractive_or_repulsive=='attractive':
                sorted_aas=sorted(use_these_aas.items(), key=lambda x: x[1])
            else:
                sorted_aas=sorted(use_these_aas.items(), key=lambda x: x[1], reverse=True)
            if len(sorted_aas) > 8:
                sorted_aas=sorted_aas[:8]
            new_sequence+=random.choice(sorted_aas)[0]

    # check to see if sequence is attractive or repulsive.
    final_eps = get_epsilon_value(new_sequence, interacting_sequence, loaded_model)
    if final_eps < 0 and attractive_or_repulsive=='repulsive':
        # try optimizing the sequence to become positive
        attempted_sequence=optimize_to_epsilon_value(new_sequence, interacting_sequence, 10,
            exclude_aas=exclude_aas, model=model, preloaded_IMC_object=loaded_model, return_best_sequence=True)
    elif final_eps > 0 and attractive_or_repulsive=='attractive':
        # try optimizing the sequence to become positive
        attempted_sequence=optimize_to_epsilon_value(new_sequence, interacting_sequence, -10,
            exclude_aas=exclude_aas, model=model, preloaded_IMC_object=loaded_model, return_best_sequence=True)
    else:
        # if we wanted it to be attractive and it was, 
        # or if we wanted it to be repulsive and it was, 
        # return it.
        return new_sequence

    # now check the attempted sequence
    attempted_eps = get_epsilon_value(attempted_sequence, interacting_sequence, loaded_model)
    if attempted_eps > 0 and attractive_or_repulsive=='repulsive':
        return attempted_sequence
    elif attempted_eps < 0 and attractive_or_repulsive=='attractive':
        return attempted_sequence
    else:
        # if we couldn't optimize it to be more attractive or repulsive,
        # raise goose_exceptions.GooseFail
        raise goose_exceptions.GooseFail(f'Failed to create {attractive_or_repulsive} sequence.')
            

def increase_epsilon(sequence_of_interest, interacting_sequence=None,
                    num_iterations=None, exclude_aas=[], model='Mpipi_GGv1',
                    maximal_optimization=False):
    '''
    Function to increase the epsilon value of a sequence relative to another sequence. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want to increase or decrease the epsilon value of.

    interacting_sequence : str
        The sequence that the sequence of interest is interacting with.
        If None, this will just do homotypic interactions (will do 
        sequence_of_interest with itself)

    num_iterations : int
        The number of iterations to run the optimization. 
        If none, default is length of sequence_of_interest/10

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    maximal_optimization : bool
        Whether to optimize to the maximal extent possible. 
        Reduces the sequence space explored, so default is False

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # load the model
    loaded_model = load_IMC_object(model)

    # if interacting sequence is None, set it to the sequence of interest
    if num_iterations is None:
        num_iterations=int(len(sequence_of_interest)/10)
    if num_iterations==0:
        num_iterations=1

    return optimize_to_epsilon_value(sequence_of_interest, interacting_sequence, 1000,
        allowed_error=0.1, optimization_iterations=num_iterations, exclude_aas=exclude_aas,
        model=model, preloaded_IMC_object=loaded_model, return_best_sequence=True,
        maximal_optimization=maximal_optimization)

def decrease_epsilon(sequence_of_interest, interacting_sequence=None,
                    num_iterations=None, exclude_aas=[], model='Mpipi_GGv1',
                    maximal_optimization=False):
    '''
    Function to decrease the epsilon value of a sequence relative to another sequence. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want to increase or decrease the epsilon value of.

    interacting_sequence : str
        The sequence that the sequence of interest is interacting with.
        If None, this will just do homotypic interactions (will do 
        sequence_of_interest with itself)

    num_iterations : int
        The number of iterations to run the optimization. 
        If none, default is length of sequence_of_interest/10

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    maximal_optimization : bool
        Whether to optimize to the maximal extent possible. 
        Reduces the sequence space explored, so default is False

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # load the model
    loaded_model = load_IMC_object(model)

    # if interacting sequence is None, set it to the sequence of interest
    if num_iterations is None:
        num_iterations=int(len(sequence_of_interest)/10)
    if num_iterations==0:
        num_iterations=1

    return optimize_to_epsilon_value(sequence_of_interest, interacting_sequence, -1000,
        allowed_error=0.1, optimization_iterations=num_iterations, exclude_aas=exclude_aas,
        model=model, preloaded_IMC_object=loaded_model, return_best_sequence=True,
        maximal_optimization=maximal_optimization)

