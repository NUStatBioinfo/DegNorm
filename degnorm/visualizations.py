from platform import platform
from re import search
import matplotlib
if search('Linux', platform()):
    matplotlib.use('agg')
import pylab as plt
import seaborn as sns