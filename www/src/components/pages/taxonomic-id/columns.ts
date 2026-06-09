import { h } from 'vue'
import type { ColumnDef, Row } from '@tanstack/vue-table'
import { ArrowUpDown, ArrowUp, ArrowDown, ChevronDown, ChevronRight } from 'lucide-vue-next'
import { Button } from '@/components/ui/button'

export interface TaxonomicIDRow {
    sample: string
    rank: number
    species: string
    ani: number
    metaSpecies: string
    metaGemsparcl: string
    metaGtdb: string
    clusterDetails?: TaxonomicIDRow[]
}

function sortableHeader(label: string, align?: 'right') {
    return ({ column }: { column: { toggleSorting: (asc: boolean) => void; getIsSorted: () => string | false } }) => {
        const sorted = column.getIsSorted()
        const icon = sorted === 'asc' ? ArrowUp : sorted === 'desc' ? ArrowDown : ArrowUpDown
        return h(Button, {
            variant: 'ghost',
            class: align === 'right' ? 'w-full justify-end' : undefined,
            onClick: () => column.toggleSorting(sorted === 'asc'),
        }, () => [label, h(icon, { class: 'ml-2 h-4 w-4' })])
    }
}

export function aniTextClass(ani: number): string {
    const percent = ani * 100
    if (percent < 75) return 'text-red-600'
    if (percent <= 85) return 'text-yellow-600'
    return ''
}

export const columns: ColumnDef<TaxonomicIDRow>[] = [
    {
        accessorKey: 'sample',
        header: sortableHeader('Sample'),
        cell: ({ row }: { row: Row<TaxonomicIDRow> }) => {
            const canExpand = row.getCanExpand()
            const isExpanded = row.getIsExpanded()
            return h('div', { class: 'flex items-center gap-1' }, [
                canExpand
                    ? h('span', {
                        class: 'h-5 w-5 p-0 text-gray-400 flex items-center justify-center',
                    }, h(isExpanded ? ChevronDown : ChevronRight, { class: 'h-3.5 w-3.5' }))
                    : h('span', { class: 'w-5' }),
                h('span', { class: 'font-medium font-mono truncate' }, row.getValue('sample')),
            ])
        },
    },
    {
        accessorKey: 'species',
        header: sortableHeader('Species'),
        cell: ({ row }) => h('div', { class: 'italic' }, row.getValue('species')),
    },
    {
        accessorKey: 'ani',
        header: sortableHeader('ANI', 'right'),
        cell: ({ row }) => {
            const ani = row.getValue('ani') as number
            return h('div', { class: ['text-right font-medium', aniTextClass(ani)] }, (ani * 100).toFixed(1) + '%')
        },
    },
    {
        accessorKey: 'metaGemsparcl',
        header: 'Gemsparcl ID',
    },
]
