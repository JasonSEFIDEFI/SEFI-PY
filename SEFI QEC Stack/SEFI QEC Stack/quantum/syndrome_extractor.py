class SyndromeExtractor:

    def extract(
        self,
        transmitted,
        received
    ):

        syndrome = []

        for i in range(len(transmitted)):

            if transmitted[i] != received[i]:
                syndrome.append(i)

        return syndrome