// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const lightCodeTheme = require('prism-react-renderer/themes/github');
const darkCodeTheme = require('prism-react-renderer/themes/dracula');

const math = require('remark-math');
const katex = require('rehype-katex');

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'Oxford Statistical Genomics Summer School',
  tagline: 'Genomics, bioinformatics and statistics training resources',
  url: 'https://well.ox.ac.uk/',
  baseUrl: '/~gav/projects/oxford_statgen_summer_school/website/',
  onBrokenLinks: 'warn',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',
  trailingSlash: true,
  // GitHub pages deployment config.
  // If you aren't using GitHub pages, you don't need these.
  organizationName: 'University of Oxford', // Usually your GitHub org/user name.
  projectName: 'statgen training', // Usually your repo name.

  // Even if you don't use internalization, you can use this field to set useful
  // metadata like html lang. For example, if your site is Chinese, you may want
  // to replace "en" with "zh-Hans".
  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          routeBasePath: '/',
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          editUrl:
            'https://github.com/gavinband/oxford-statgen-summer-school/edit/master',
          remarkPlugins: [math],
          rehypePlugins: [katex],
        },
        blog: false,
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      }),
    ],
  ],
  stylesheets: [
    {
      href: 'https://cdn.jsdelivr.net/npm/katex@0.13.24/dist/katex.min.css',
      type: 'text/css',
      integrity:
        'sha384-odtC+0UGzzFL/6PNoE8rX/SPcQDXBJ+uRepguP4QkPCm2LBxH3FA3y+fKSiJ+AmM',
      crossorigin: 'anonymous',
    },
  ],
  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: 'Oxford Statistical Genomics Summer School',
        logo: {
          alt: 'Statistical Genomics Summer School',
          src: 'img/statistics-summer-school-logo.jpg',
        },
        items: [
          {
            type: 'doc',
            docId: 'Overview',
            position: 'left',
            label: 'Practicals',
          },
          {
            href: 'https://github.com/gavinband/oxford-statgen-summer-school',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
	  footer: {
        style: 'dark',
        links: [
          {
            title: 'Docs',
            items: [
                {
                    label: 'Contributors',
                    href: 'contributors'
                }
            ],
          },
          {
            title: 'More',
            items: [
              {
                label: 'GitHub',
                href: 'https://github.com/gavinband/oxford-statgen-summer-school',
              },
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} University of Oxford.  Built with Docusaurus.`,
      },
      prism: {
        theme: lightCodeTheme,
        darkTheme: darkCodeTheme,
      },
    }),
};

module.exports = config;
