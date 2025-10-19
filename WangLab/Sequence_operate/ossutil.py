import argparse
import os


def prepare_argparser(argparser):
    argparser.add_argument("-conf", '--config', dest='config', type=str, required=True, nargs="?",
                           help='Config file')

    return


def download(args):
    with open(args.config) as fi:
        while True:
            content = fi.readline()
            if not content:
                break
            if content.startswith("#"):
                continue
            # check
            paras = content.rstrip('\n').split('\t')
            if not len(paras) == 5:
                print(f'Line "{content}" missing some parameters, please check.')
                continue

            keyid, key, path, e, out_dir = paras
            cmd = f'{os.path.join("/home/byf1999/Softwares", "ossutil64")} cp -e {e} -i {keyid} -k {key} -r {path} {out_dir}'
            os.system(cmd)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    prepare_argparser(argparser)
    args = argparser.parse_args()
    download(args)
