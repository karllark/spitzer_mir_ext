import argparse

from measure_extinction.stardata import StarData
from measure_extinction.extdata import ExtData


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("redstarname", help="name of reddened star")
    parser.add_argument("compstarname", help="name of comparision star")
    parser.add_argument(
        "--path",
        help="base path to observed data",
        default="/home/kgordon/Python_git/extstar_data/",
    )
    args = parser.parse_args()

    # read in the observed data for both stars
    redstarobs = StarData("DAT_files/%s.dat" % args.redstarname, path=args.path)
    compstarobs = StarData("DAT_files/%s.dat" % args.compstarname, path=args.path)

    # output filebase
    filebase = "fits/%s_%s" % (args.redstarname, args.compstarname)

    # calculate the extinction curve
    extdata = ExtData()
    extdata.calc_elx(redstarobs, compstarobs)

    # save the extinction curve
    out_fname = "fits/%s_%s_ext.fits" % (args.redstarname, args.compstarname)
    extdata.save(out_fname)
