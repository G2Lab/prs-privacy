ALL_FIELDS = ["rsID",
              "chr_name",
              "chr_position",
              "effect_allele",
              "other_allele",
              "effect_weight",
              "allelefrequency_effect",
              "locus_name",
              "variant_description",
              "OR",
              "hm_source",
              "hm_rsID",
              "hm_chr",
              "hm_pos",
              "hm_inferOtherAllele"]


class Variant:
    def __init__(self, fields=None):
        if (fields is not None) and (fields['effect_weight'] is not None):
            fields['effect_weight'] = float(fields['effect_weight'])
        for field in ALL_FIELDS:
            if field not in fields:
                fields[field] = None
        self.fields = fields

    def get_hm_chr(self):
        return self.fields['hm_chr']

    def get_hm_pos(self):
        return self.fields['hm_pos']

    def get_id(self):
        return self.fields['rsID']

    def get_weight(self):
        return self.fields['effect_weight']


class PGS:
    def __init__(self, pgs_id=None, trait_name=None, trait_efo=None,
                 genome_build=None, weight_type=None, HmPOS_build=None,
                 variants_number=None, fieldnames=None):
        self.pgs_id = pgs_id
        self.trait_name = trait_name
        self.trait_efo = trait_efo
        self.genome_build = genome_build
        self.weight_type = weight_type
        self.HmPOS_build = HmPOS_build
        self.variants_number = variants_number
        self.fieldnames = fieldnames
        self.variants = {}
        self.chr_to_weight = {}

    def load(self, input_file):
        with open(input_file, "r") as f:
            for line in f:
                line = line.rstrip()
                #  If it is the header
                if line.startswith("#"):
                    if line[1:].startswith("pgs_id"):
                        self.pgs_id = line.split("=")[1]
                    elif line[1:].startswith("trait_mapped"):
                        self.trait_name = line.split("=")[1]
                    elif line.startswith("trait_efo"):
                        self.trait_efo = line.split("=")[1]
                    elif line[1:].startswith("genome_build"):
                        self.genome_build = line.split("=")[1]
                    elif line[1:].startswith("weight_type"):
                        self.weight_type = line.split("=")[1]
                    elif line[1:].startswith("HmPOS_build"):
                        self.HmPOS_build = line.split("=")[1]
                    elif line[1:].startswith("variants_number"):
                        self.variants_number = int(line.split("=")[1])
                    continue

                # Specified variant fields
                if line.startswith("rsID"):
                    self.fieldnames = line.split("\t")
                    continue

                fields = dict(zip(self.fieldnames, line.split("\t")))
                variant = Variant(fields)
                self.variants[variant.get_id()] = variant
                self.chr_to_weight[f'{variant.get_hm_chr()}:{variant.get_hm_pos()}'] = variant.get_weight()
