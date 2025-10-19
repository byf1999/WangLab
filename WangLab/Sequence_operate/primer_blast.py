"""
This script read sequence file and automate input into NCBI primer-blast program
"""
import time

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

from Bio import SeqIO
import argparse
import sys
import re
import os


def prepare_argparser(arg_parser):
    arg_parser.add_argument("-i", '--input', dest='input', type=str, required=True, nargs="?", help='Input sequences')
    arg_parser.add_argument("-r", '--ref', dest='ref', type=str, required=False, nargs="?", help='Reference sequence',
                            default=None)
    arg_parser.add_argument("-P", '--p5', dest='p5', type=str, required=False, nargs=2, help="5' primer range",
                            default=None)
    arg_parser.add_argument("-p", '--p3', dest='p3', type=str, required=False, nargs=2, help="3' primer range",
                            default=None)
    arg_parser.add_argument("-s", '--size', dest='product_size', type=str, required=False, nargs=2,
                            help='Product size range, default is 70-1000 (same with NCBI default)',
                            default=['70', '1000'])
    arg_parser.add_argument("-t", '--tm', dest='TM', type=str, required=False, nargs=3,
                            help='Primer Tm range, min, opt, and max', default=['60', '62', '63'])
    arg_parser.add_argument('--tm_diff', dest='PRIMER_MAX_DIFF_TM', type=str, required=False, nargs="?",
                            help='PRIMER_MAX_DIFF_TM', default='3')
    arg_parser.add_argument('--explorer', dest='explorer', type=str, required=False, nargs="?",
                            help='Explorer, default is Firefox', default='Firefox')
    arg_parser.add_argument('--num_primer', dest='PRIMER_NUM_RETURN', type=str, required=False, nargs="?",
                            help='PRIMER_NUM_RETURN', default='10')
    arg_parser.add_argument("-o", '--output', dest='output', type=str, required=False, nargs="?",
                            help='Output file path, default is stdout', default=None)
    arg_parser.add_argument("-w", '--window', dest='window', type=str, required=False, nargs="?", default=False,
                            help='Whether close blast window, default is False, in this mode, program will pause after each primer pair is designed.')
    # arg_parser.add_argument("-i", '--input', dest='input', type=str, required=True, nargs="?", help='Input sequences')


def primer_blast_single_sequence(input_sequence, CUSTOMSEQFILE, PRIMER5, PRIMER3, PRIMER_PRODUCT_MIN,
                                 PRIMER_PRODUCT_MAX, PRIMER_MIN_TM,
                                 PRIMER_OPT_TM, PRIMER_MAX_TM, PRIMER_MAX_DIFF_TM=3, PRIMER_NUM_RETURN=10,
                                 explorer='Firefox'):
    # 创建浏览器对象
    EXPLORER = {'Firefox': webdriver.Firefox,
                'Chrome': webdriver.Chrome,
                'Edge': webdriver.ChromiumEdge}
    driver = EXPLORER[explorer]()

    # 打开网页
    driver.get("https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=BlastHome")

    # seq = driver.find_element(By.NAME, "INPUT_SEQUENCE")
    # seq.clear()
    # seq.send_keys(INPUT_SEQUENCE)

    # 用JavaScript来设置输入框的value属性为你想要输入的序列
    driver.execute_script(f"document.getElementsByName('INPUT_SEQUENCE')[0].value='{input_sequence}'")

    if PRIMER5:
        PRIMER5_START, PRIMER5_END = PRIMER5
        primer5_start = driver.find_element(By.NAME, "PRIMER5_START")
        primer5_start.clear()
        primer5_start.send_keys(PRIMER5_START)

        primer5_end = driver.find_element(By.NAME, "PRIMER5_END")
        primer5_end.clear()
        primer5_end.send_keys(PRIMER5_END)
    if PRIMER3:
        PRIMER3_START, PRIMER3_END = PRIMER3

        if eval(PRIMER3_START) <= 0:
            PRIMER3_START = str(len(input_sequence) + eval(PRIMER3_START))
        if eval(PRIMER3_END) <= 0:
            PRIMER3_END = str(len(input_sequence) + eval(PRIMER3_END))

        primer3_start = driver.find_element(By.NAME, "PRIMER3_START")
        primer3_start.clear()
        primer3_start.send_keys(PRIMER3_START)

        primer3_end = driver.find_element(By.NAME, "PRIMER3_END")
        primer3_end.clear()
        primer3_end.send_keys(PRIMER3_END)

    primer_product_min = driver.find_element(By.NAME, "PRIMER_PRODUCT_MIN")
    primer_product_min.clear()
    primer_product_min.send_keys(PRIMER_PRODUCT_MIN)

    primer_product_max = driver.find_element(By.NAME, "PRIMER_PRODUCT_MAX")
    primer_product_max.clear()
    primer_product_max.send_keys(PRIMER_PRODUCT_MAX)

    primer_num_return = driver.find_element(By.NAME, "PRIMER_NUM_RETURN")
    primer_num_return.clear()
    primer_num_return.send_keys(PRIMER_NUM_RETURN)

    primer_min_tm = driver.find_element(By.NAME, "PRIMER_MIN_TM")
    primer_min_tm.clear()
    primer_min_tm.send_keys(PRIMER_MIN_TM)

    primer_opt_tm = driver.find_element(By.NAME, "PRIMER_OPT_TM")
    primer_opt_tm.clear()
    primer_opt_tm.send_keys(PRIMER_OPT_TM)

    primer_max_tm = driver.find_element(By.NAME, "PRIMER_MAX_TM")
    primer_max_tm.clear()
    primer_max_tm.send_keys(PRIMER_MAX_TM)

    primer_max_diff_tm = driver.find_element(By.NAME, "PRIMER_MAX_DIFF_TM")
    primer_max_diff_tm.clear()
    primer_max_diff_tm.send_keys(PRIMER_MAX_DIFF_TM)

    wait = WebDriverWait(driver, 10)
    wait.until(EC.presence_of_element_located((By.NAME, "PRIMER_SPECIFICITY_DATABASE")))
    PRIMER_SPECIFICITY_DATABASE = Select(driver.find_element(By.NAME, "PRIMER_SPECIFICITY_DATABASE"))
    PRIMER_SPECIFICITY_DATABASE.select_by_value("Custom")

    upload = driver.find_element(By.NAME, "CUSTOMSEQFILE")
    upload.send_keys(CUSTOMSEQFILE)

    # 等待按钮可点击
    btnSubmit = wait.until(EC.element_to_be_clickable((By.CSS_SELECTOR, "[value='Get Primers']")))
    # 点击按钮
    btnSubmit.click()
    return driver


def select_primer(driver):
    page_source = driver.page_source
    # r = re.search('<div class="prPairInfo">(.*)<!--/#alignments-->', page_source, flags=re.S)
    # r = re.search(r'<div class="prPairInfo">(.*)<div class="prPairDtl">', page_source, flags=re.S)
    r = re.findall(r'<div class="prPairInfo">(.*?)<div class="prPairDtl">', page_source, flags=re.S)
    while not r:
        time.sleep(5)  # try load primer data every 5 seconds
        page_source = driver.page_source
        r = re.findall(r'<div class="prPairInfo">(.*?)<div class="prPairDtl">', page_source, flags=re.S)

    best_s = 40
    best_primers = []
    for i, align in enumerate(r):
        forward = re.search('<th>Forward primer</th>(.*?)</tr>', align).group(1)
        tmp_f = re.findall('<td>(.*?)</td>', forward)
        reverse = re.search('<th>Reverse primer</th>(.*?)</tr>', align).group(1)
        tmp_r = re.findall('<td>(.*?)</td>', reverse)
        num = align.count('Product length')
        # len(re.findall('Product length', align))
        s = sum([eval(i) for i in tmp_f[7:] + tmp_r[7:]])
        if num == 1 and best_s > s:
            best_s = s
            best_primers = [tmp_f[0], tmp_r[0]]

    return best_primers


def primer_blast(args):
    # print(args)
    input_file = args.input
    ref_file = args.ref
    PRIMER5 = args.p5
    PRIMER3 = args.p3
    PRIMER_PRODUCT_MIN, PRIMER_PRODUCT_MAX = args.product_size
    PRIMER_MIN_TM, PRIMER_OPT_TM, PRIMER_MAX_TM = args.TM
    PRIMER_MAX_DIFF_TM = args.PRIMER_MAX_DIFF_TM
    PRIMER_NUM_RETURN = args.PRIMER_NUM_RETURN
    explorer = args.explorer
    file = args.output
    window = args.window
    if file and os.path.exists(file):
        print(f"{file} already exists!")
        return

    primers = {}
    for record in SeqIO.parse(input_file, 'fasta'):
        print(f'Dealing with {record.id}:')
        # print(record.seq, ref_file, PRIMER5, PRIMER3, PRIMER_PRODUCT_MIN, PRIMER_PRODUCT_MAX, PRIMER_MIN_TM,
        #       PRIMER_OPT_TM, PRIMER_MAX_TM, PRIMER_MAX_DIFF_TM, PRIMER_NUM_RETURN, explorer)
        driver = primer_blast_single_sequence(record.seq, ref_file, PRIMER5, PRIMER3, PRIMER_PRODUCT_MIN,
                                              PRIMER_PRODUCT_MAX, PRIMER_MIN_TM,
                                              PRIMER_OPT_TM, PRIMER_MAX_TM, PRIMER_MAX_DIFF_TM=PRIMER_MAX_DIFF_TM,
                                              PRIMER_NUM_RETURN=PRIMER_NUM_RETURN, explorer=explorer)
        primer = select_primer(driver)
        primers[record.id] = primer
        if window:
            driver.close()
        else:
            input('Press any key to continue:')

    if file:
        with open(file, 'w') as fo:
            for k, v in primers.items():
                fo.write('\t'.join([k, *v]) + '\n')
    else:
        for k, v in primers.items():
            print('\t'.join([k, *v]))


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    prepare_argparser(arg_parser)
    args = arg_parser.parse_args()

    primer_blast(args)
