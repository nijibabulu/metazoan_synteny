#! /usr/bin/env python

import argparse

import openpyxl


def takewhile_inclusive(pred, i):
    # takewhile_inclusive(lambda x: x<5, [1,4,6,4,1]) --> 1 4 5
    for x in i:
        yield x
        if not pred(x):
            break


class CellIdentity(object):
    def __init__(self, last_worksheet, name=None):
        self.last_worksheet = last_worksheet
        self.name = name or 'Unnamed {}'.format(last_worksheet)
        self.metacells = []

    def parse(self, ws):
        our_worksheets = takewhile_inclusive(
            lambda w: w.title != 'C{}'.format(self.last_worksheet), ws)
        self.metacells.extend(w.title[1:] for w in our_worksheets)


def handle_identities(worksheet_breaks, xlsx_filename, out_filename):
    # type: (list[tuple], str, str) -> None
    cell_identities = map(lambda t: CellIdentity(*t), worksheet_breaks)
    xlsx = openpyxl.load_workbook(xlsx_filename)

    ws_iter = iter(xlsx.worksheets)
    map(lambda ci: ci.parse(ws_iter), cell_identities)

    cell_identity_map = {mc: c.name
                         for c in cell_identities
                         for mc in c.metacells}

    with open(out_filename, 'w') as f:
        f.write('\n'.join(
            '{}\t{}'.format(mc, ci) for mc, ci in cell_identity_map.items()))


if __name__ == '__main__':
    """
    The metacell identities are not organized in order of their cell identity 
    clusters, but by their ordering in the excel sheet. In order to avoid 
    manual copying, we created this script which goes through the worksheets
    in order. 
    """

    aq_adult_breaks = [('4', 'Choanocytes 1'),
                       ('14', 'Choanocytes 2'),
                       ('23',),
                       ('45', 'Pinacocytes 1'),
                       ('52', 'Pinacocytes 2'),
                       ('34', 'Archaeocytes 1'),
                       ('28', 'Archaeocytes 2'),
                       ('27',),
                       ('26', 'Aspzincin Cells'),
                       ('22', 'Collagen Cells'),
                       ('21', 'Sperm'),
                       ('25', 'Host Defense Cells'),
                       ('24', 'Unnamed 24')]

    aq_larva_breaks = [('3', 'Ciliated epithelium'),
                        ('1', 'Posterior pole cells'),
                        ('2', 'Flask cells'),
                        ('15', 'Non-ciliated epithlium'),
                        ('13',), ('14',),
                        ('21', 'Archaeocyte-like cells'),
                        ('12', 'Anterior pole cells')]

    ml_breaks = [('10', 'Digestive cells'),
                 ('17', 'Epithelial cells'),
                 ('27',), ('28',), ('26',), ('44',),
                 ('43', 'Smooth muscle'),
                 ('47', 'Striated muscle'),
                 ('50', 'Comb cells'),
                 ('42',), ('39',), ('41',), ('38',), ('40',), ('32',), ('33',),
                 ('34',), ('55',), ('35',), ('54',),
                 ('52', 'Shk protein-producing cells'),
                 ('30',), ('29',), ('53',),
                 ('51', 'Photocytes',),
                 ('37',), ('31',), ('36',)]

    ta_breaks = [('11', 'Lipiophil cells'),
                 ('1',),
                 ('22', 'Digestive cells'),
                 ('31', 'Fibre cells'),
                 ('37',),
                 ('43', 'Epithelial cells'),
                 ('36', 'Peptidergic cells'),
                 ('39',)]

    parser = argparse.ArgumentParser(
        description='Extract the metacell ids as they are in Figure 1c of '
                    'doi: 10.1038/s41559-018-0575-6')
    parser.add_argument('XLSX_AQ_ADULT')
    parser.add_argument('XLSX_AQ_LARVA')
    parser.add_argument('XLSX_ML')
    parser.add_argument('XLSX_TA')
    parser.add_argument('AQ_ADULT_OUT')
    parser.add_argument('AQ_LARVA_OUT')
    parser.add_argument('ML_OUT')
    parser.add_argument('TA_OUT')
    args = parser.parse_args()

    handle_identities(aq_adult_breaks, args.XLSX_AQ_ADULT, args.AQ_ADULT_OUT)
    handle_identities(aq_larva_breaks, args.XLSX_AQ_LARVA, args.AQ_LARVA_OUT)
    handle_identities(ml_breaks, args.XLSX_ML, args.ML_OUT)
    handle_identities(ta_breaks, args.XLSX_TA, args.TA_OUT)

