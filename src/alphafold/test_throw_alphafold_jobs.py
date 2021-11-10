import pytest

from throw_alphafold_jobs import DetermineJobParams

@pytest.mark.parametrize(('length', 'estimated_time', 'ensemble'), [
    (100, '3:29:00', 'full'),
    (116, '4:01:00', 'full'),
    (145, '5:00:00', 'full'),
    (521, '3:31:00', 'noens')
])
def test_estimated_runtime(length, estimated_time, ensemble):
    assert DetermineJobParams.estimate_time_from_length(length, ensemble) == estimated_time


def main():
    samples = [
        (100, 'full'),
        (200, 'full'),
        (300, 'full'),
        (350, 'full'),
        (400, 'noens'),
        (500, 'noens'),
        (600, 'noens'),
        (700, 'noens')
    ]
    for sample in samples:
        l, e = sample
        print(l, e, DetermineJobParams.estimate_time_from_length(l, e))


if __name__ == '__main__':
    main()
