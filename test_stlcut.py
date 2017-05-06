import admesh
import glob
import math
import os
import pytest
import subprocess
import random

OUT1 = 'Cut_Mesh_1.stl'
OUT2 = 'Cut_Mesh_2.stl'

STLS = glob.glob('stl_files/*.stl')

def rm_outputs():
    for out in [OUT1, OUT2]:
        try:
            os.remove(out)
        except FileNotFoundError:
            pass


def cut(path, *args):
    rm_outputs()
    args = [str(a) for a in args]
    subprocess.call(['./cut2', path] + args)


def volume(path):
    stl = admesh.Stl(path)
    stl.calculate_volume()
    return stl.stats['volume']


def volume_parts(part1=OUT1, part2=OUT2):
    try:
        volume_1 = volume(part1)
    except admesh.AdmeshError:
        volume_1 = 0
    try:
        volume_2 = volume(part2)
    except admesh.AdmeshError:
        volume_2 = 0
    return volume_1 + volume_2


@pytest.mark.parametrize('stl', STLS)
def test_cut_z0_sum_of_volumes_match_original_volume(stl):
    if stl == 'stl_files/duplo.stl':
        pytest.xfail('Known issue with p2t')
    cut(stl, 0, 0, 1, 0.5)
    volume_total = volume(stl)
    volume_sum = volume_parts()
    assert math.isclose(volume_sum, volume_total, rel_tol=0.05)
